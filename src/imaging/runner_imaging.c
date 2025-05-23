/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Will Roper (w.roper@sussex.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include <config.h>

/* Standard includes. */
#include <math.h>

/* This files header. */
#include "images.h"

/* Local includes. */
#include "cell.h"
#include "lightcone/projected_kernel.h"
#include "runner.h"
#include "timers.h"

#include <math.h>

/**
 * @brief Loop over all pixels within a particle's kernel support radius
 * and add the kernel value to the image.
 *
 * This is a generic function for smoothing any quantity into an image based
 * on a particle's position and a kernel.
 *
 * @param image  The image to be smoothed into.
 * @param x     The particle's position in world coordinates.
 * @param nx    The number of pixels in the X direction.
 * @param ny    The number of pixels in the Y direction.
 * @param dx    The pixel size in the X direction.
 * @param dy    The pixel size in the Y direction.
 * @param x0    The lower-left corner of the image in world coordinates.
 * @param y0    The lower-left corner of the image in world coordinates.
 * @param h     The kernel width in world units.
 * @param value The value to be added to the image.
 * @param kt    The kernel type.
 * @param p     The particle to be smoothed.
 */
static void runner_do_kernel_smoothed_image_loop(
    double *image, double x[3], int nx, int ny, double dx, double dy, double x0,
    double y0, double h, double value, struct projected_kernel_table *kt) {

  /* Particle position in pixel coords */
  double px = (x[0] - x0) / dx;
  double py = (x[1] - y0) / dy;

  /* Particle position in pixel coords */
  int ix0 = (int)floor(px);
  int iy0 = (int)floor(py);

  /* Kernel support radius: u_max * h (in world units) */
  double support = kt->u_max * h;
  int delta_pix = (int)ceil(support / dx);

  /* Loop over all pixels within ±delta_pix of the particle. */
  for (int di = -delta_pix; di <= delta_pix; di++) {
    int ix = ix0 + di;

    /* Skip out-of-bounds pixels. */
    if (ix < 0 || ix >= nx)
      continue;

    /* Compute the distance from the pixel center to the particle. */
    double rx = fabs(px - (ix + 0.5)) * dx;

    /* Loop over all pixels in the Y direction. */
    for (int dj = -delta_pix; dj <= delta_pix; dj++) {
      int iy = iy0 + dj;

      /* Skip out-of-bounds pixels. */
      if (iy < 0 || iy >= ny)
        continue;

      /* Compute the distance from the pixel center to the particle. */
      double ry = fabs(py - (iy + 0.5)) * dy;

      /* Compute the distance from the pixel center to the particle. */
      double b = sqrt(rx * rx + ry * ry);

      /* Compute the impact parameter. */
      double u = b / h;

      /* Skip out-of-bounds pixels. */
      if (u >= kt->u_max)
        continue;

      /* Compute the kernel value. */
      double w = projected_kernel_eval(kt, u);
      image[ix * ny + iy] += value * w;
    }
  }
}

static double runner_get_part_weight(struct image_data *image_data,
                                     struct part *p) {

  /* Get the weight for this particle. */
  if (strcmp(image_data->field_name, "mass") == 0) {
    return p->mass;
  } else if (strcmp(image_data->field_name, "density") == 0) {
    return p->density.wcount + 32.0 / 3;
  } else if (strcmp(image_data->field_name, "temperature") == 0) {
    return p->cooling_data.subgrid_temp;
  } else {
    error("Unknown field name for gas: %s", image_data->field_name);
  }
}

/**
 * @brief Make a smoothed mass image for dark matter.
 *
 * Each dark-matter particle’s mass is distributed over nearby pixels
 * using a 2D Gaussian of width ε (the softening length), truncated
 * at 3 ε for efficiency.
 *
 * @param image_data The overall image settings (pixel size, etc.).
 * @param iimage     Index of which image buffer to write into.
 * @param c          Cell containing the dark-matter particles.
 */
static void
runner_do_gpart_mass_image(struct image_common_data *image_common_data,
                           int iimage, struct cell *c) {
  /* Extract the image data. */
  struct image_data *image_data = &image_common_data->images[iimage];

  /* Extract the particle list and count from the cell’s gravity data. */
  struct gpart *gparts = c->grav.parts;
  int gcount = c->grav.count;

  /* Origin (lower-left) of this cell in “world” coordinates. */
  double x0 = c->image_data.padded_loc[0];
  double y0 = c->image_data.padded_loc[1];

  /* Pixel dimensions in world units. */
  double dx = image_data->pixel_size[0];
  double dy = image_data->pixel_size[1];

  /* Number of pixels along X and Y within this padded cell. */
  int nx = (int)ceil(c->image_data.padded_width[0] / dx);
  int ny = (int)ceil(c->image_data.padded_width[1] / dy);

  /* Pointer to the raw image buffer we’ll accumulate into. */
  double *image = c->image_data.images[iimage];

  /*
   * Pre-declare a small buffer for weights. We assume
   * 3ε/dx will never exceed 15 pixels or so; adjust size
   * if you use a larger cutoff.
   */
  const int MAX_RAD_PX = 20;
  double wbuf[2 * MAX_RAD_PX + 1][2 * MAX_RAD_PX + 1];

  /* Loop over each dark-matter “particle.” */
  for (int i = 0; i < gcount; i++) {
    struct gpart *gp = &gparts[i];
    double eps = gp->epsilon; /* softening length */

    /*
     * Convert particle position into floating‐point pixel coords:
     *   fx = (world_x – x0) / dx
     *   fy = (world_y – y0) / dy
     */
    double fx = (gp->x[0] - x0) / dx;
    double fy = (gp->x[1] - y0) / dy;

    /*
     * Compute cutoff radius in pixel units (3 × ε / dx).
     * We’ll loop over all pixels within ±delta_pix of the particle.
     */
    double eps_px = eps / dx;
    double radius = 3.0 * eps_px;
    int delta_pix = (int)ceil(radius);

    /*
     * Sanity: clamp to our buffer size.
     * If you ever need a bigger radius, enlarge MAX_RAD_PX.
     */
    if (delta_pix > MAX_RAD_PX) {
      delta_pix = MAX_RAD_PX;
    }

    /*
     * 1) Compute un-normalized Gaussian weights over the patch.
     *    W(r) = (1 / (2π ε²)) exp[–½ (r/ε)²], with r² = dx²+dy².
     *    Here we measure distances in pixel units.
     */
    double total_w = 0.0;
    double ix0 = floor(fx);
    double iy0 = floor(fy);
    for (int di = -delta_pix; di <= delta_pix; di++) {
      for (int dj = -delta_pix; dj <= delta_pix; dj++) {
        /*
         * Center each pixel at (ix0 + 0.5 + di, iy0 + 0.5 + dj):
         * compute offset from true particle position.
         */
        double rx = fx - (ix0 + 0.5 + di);
        double ry = fy - (iy0 + 0.5 + dj);
        double r2 = rx * rx + ry * ry;

        /* Gaussian in pixel units: σ = eps_px */
        double w =
            exp(-0.5 * r2 / (eps_px * eps_px)) / (2.0 * M_PI * eps_px * eps_px);

        wbuf[di + delta_pix][dj + delta_pix] = w;
        total_w += w;
      }
    }

    /*
     * 2) Normalize and distribute mass:
     *    each pixel gets gp->mass * (w / total_w)
     */
    int base_ix = (int)ix0 - delta_pix;
    int base_iy = (int)iy0 - delta_pix;
    for (int di = 0; di <= 2 * delta_pix; di++) {
      int ix = base_ix + di;
      if (ix < 0 || ix >= nx)
        continue; /* skip out-of-bounds */
      for (int dj = 0; dj <= 2 * delta_pix; dj++) {
        int iy = base_iy + dj;
        if (iy < 0 || iy >= ny)
          continue;

        double w_norm = wbuf[di][dj] / total_w;
        image[ix * ny + iy] += gp->mass * w_norm;
      }
    }
  }
}

/**
 * @brief Mass-weighted image for gas particles using SPH smoothing.
 *
 * Each part is smoothed into the image over its SPH kernel.
 *
 * @param image_data The overall image settings (pixel size, etc.).
 * @param iimage     Index of which image buffer to write into.
 * @param c          Cell containing the gas particles.
 */
static void
runner_do_part_mass_image(struct image_common_data *image_common_data,
                          int iimage, struct cell *c) {

  /* Extract the image data. */
  struct image_data *image_data = &image_common_data->images[iimage];

  /* Extract the particles. */
  struct part *parts = c->hydro.parts;
  int gcount = c->hydro.count;

  /* Get the image data for convenience. */
  double *image = c->image_data.images[iimage];
  double x0 = c->image_data.padded_loc[0];
  double y0 = c->image_data.padded_loc[1];
  double dx = image_data->pixel_size[0];
  double dy = image_data->pixel_size[1];
  int nx = c->image_data.num_pixels[0];
  int ny = c->image_data.num_pixels[1];

  /* Get the kernel look up table. */
  struct projected_kernel_table *kt = image_common_data->projected_kernel_table;

  /* Loop over each particle. */
  for (int i = 0; i < gcount; i++) {

    /* Get the particle and its smoothing length and mass. */
    struct part *p = &parts[i];
    double h = p->h;
    double mass = p->mass;

    /* Adding this particles contribution to the image. */
    runner_do_kernel_smoothed_image_loop(image, p->x, nx, ny, dx, dy, x0, y0, h,
                                         mass, kt);
  }
}

/**
 * @brief Temperature image for gas particles using SPH smoothing.
 *
 * Each part is smoothed into the image over its SPH kernel.
 *
 * @param image_data The overall image settings (pixel size, etc.).
 * @param iimage     Index of which image buffer to write into.
 * @param c          Cell containing the gas particles.
 */
static void
runner_do_part_temperature_image(struct image_common_data *image_common_data,
                                 int iimage, struct cell *c) {

  /* Extract the image data. */
  struct image_data *image_data = &image_common_data->images[iimage];

  /* Extract the particles. */
  struct part *parts = c->hydro.parts;
  int gcount = c->hydro.count;

  /* Get the image data for convenience. */
  double *image = c->image_data.images[iimage];
  double x0 = c->image_data.padded_loc[0];
  double y0 = c->image_data.padded_loc[1];
  double dx = image_data->pixel_size[0];
  double dy = image_data->pixel_size[1];
  int nx = c->image_data.num_pixels[0];
  int ny = c->image_data.num_pixels[1];

  /* Get the kernel look up table. */
  struct projected_kernel_table *kt = image_common_data->projected_kernel_table;

  /* Loop over each particle. */
  for (int i = 0; i < gcount; i++) {

    /* Get the particle and its smoothing length and temperature. */
    struct part *p = &parts[i];
    double h = p->h;
    double temp = p->cooling_data.subgrid_temp;

    /* Adding this particles contribution to the image. */
    runner_do_kernel_smoothed_image_loop(image, p->x, nx, ny, dx, dy, x0, y0, h,
                                         temp, kt);
  }
}

/**
 * @brief Mass-weighted temperature image for gas particles using SPH smoothing.
 *
 * Each part is smoothed into the image over its SPH kernel.
 *
 * The temperatures will be weigthed by the particle mass. This means we will
 * make two images, one for the mass and one for the temperature * mass, before
 * dividing out the mass dependence.
 *
 * @param image_common_data The overall image settings (pixel size, etc.).
 * @param iimage     Index of which image buffer to write into.
 * @param c          Cell containing the gas particles.
 */
static void runner_do_part_weighted_temperature_image(
    struct image_common_data *image_common_data, int iimage, struct cell *c) {

  /* Extract the image data. */
  struct image_data *image_data = &image_common_data->images[iimage];

  /* Extract the weight image information. */
  struct image_data *weight_image_data =
      &image_common_data->images[image_data->weight_by];

  /* Extract the particles. */
  struct part *parts = c->hydro.parts;
  int gcount = c->hydro.count;

  /* Get the image data for convenience. */
  double *image = c->image_data.images[iimage];
  double x0 = c->image_data.padded_loc[0];
  double y0 = c->image_data.padded_loc[1];
  double dx = image_data->pixel_size[0];
  double dy = image_data->pixel_size[1];
  int nx = c->image_data.num_pixels[0];
  int ny = c->image_data.num_pixels[1];

  /* Get the kernel look up table. */
  struct projected_kernel_table *kt = image_common_data->projected_kernel_table;

  /* Loop over each particle. */
  for (int i = 0; i < gcount; i++) {

    /* Get the particle and its smoothing length and temperature. */
    struct part *p = &parts[i];
    double h = p->h;
    double temp = p->cooling_data.subgrid_temp;
    double weight = runner_get_part_weight(weight_image_data, p);

    /* Adding this particles contribution to the image. */
    runner_do_kernel_smoothed_image_loop(image, p->x, nx, ny, dx, dy, x0, y0, h,
                                         temp * weight, kt);
  }
}

/**
 * @brief Compute the images for each particle in this cell.
 *
 * @param r The runner data structure.
 * @param c The cell data structure.
 * @param timer Are we timing this?
 */
void runner_do_imaging(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  /* Get the image data for convenience. */
  struct image_common_data *image_data = r->e->image_data;
  int nimages = image_data->num_images;

  /* Do we need to init the images? */
  if (c->image_data.images == NULL && cell_init_images(c, r->e) != 0) {
    error("Failed to allocate memory for the images.");
    return;
  }

  /* Zero this cells images. */
  for (int i = 0; i < nimages; i++) {
    bzero(c->image_data.images[i], c->image_data.num_pixels[0] *
                                       c->image_data.num_pixels[1] *
                                       sizeof(double));
  }

  /* Loop over all images. */
  for (int i = 0; i < nimages; i++) {

    /* What particle type are we working on? */
    const int ptype = image_data->images[i].particle_type;

    /* What field are we making an image of? */
    const char *field_name = image_data->images[i].field_name;

    /* Are we weighting the image by another image? */
    const int weight_by = image_data->images[i].weight_by;

    /* Choose the right function for the job. */
    switch (ptype) {
    case 0: /* Gas */
      if (strcmp(field_name, "mass") == 0) {
        runner_do_part_mass_image(image_data, i, c);
        // } else if (strcmp(field_name, "density") == 0) {
        //   runner_do_part_density_image(&image_data, i, c);
      } else if (strcmp(field_name, "temperature") == 0 && weight_by == -1) {
        runner_do_part_temperature_image(image_data, i, c);
      } else if (strcmp(field_name, "temperature") == 0 && weight_by != -1) {
        runner_do_part_weighted_temperature_image(image_data, i, c);
      } else {
        error("Unknown field name for gas: %s", field_name);
      }
      break;
    case 1: /* Dark matter */
      if (strcmp(field_name, "mass") == 0) {
        runner_do_gpart_mass_image(image_data, i, c);
      } else {
        error("Unknown field name for dark matter: %s", field_name);
      }
      break;
    case 2: /* Nothing */
      error("PartType 2 is not supported for imaging.");
      break;
    case 3: /* Nothing */
      error("PartType 3 is not supported for imaging.");
      break;
    case 4: /* Stars */
      error("Stars imaging not implemented yet.");
      // if (strcmp(field_name, "mass") == 0) {
      //   runner_do_spart_mass_image(&image_data, i, c);
      // } else {
      //   error("Unknown field name for stars: %s", field_name);
      // }
      break;
    case 5: /* Black holes */
      error("Black holes imaging not implemented yet.");
      // if (strcmp(field_name, "mass") == 0) {
      //   runner_do_bpart_mass_image(&image_data, i, c);
      // } else {
      //   error("Unknown field name for black holes: %s", field_name);
      // }
      break;
    default:
      error("Unknown particle type: %d", ptype);
    }
  }

  if (timer)
    TIMER_TOC(timer_dograv_down);
}

static void runner_do_recursive_imaging_collect(struct runner *r,
                                                struct cell *c,
                                                double *top_image,
                                                const int iimage) {

  /* Recurse into the progeny */
  for (int i = 0; i < 8; i++) {
    if (c->progeny[i] != NULL) {
      runner_do_recursive_imaging_collect(r, c->progeny[i], top_image, iimage);
    }
  }

  /* Get the image data for convenience. */
  const struct image_common_data *image_data = r->e->image_data;
  const double pixel_size[2] = {image_data->pixel_size[0],
                                image_data->pixel_size[1]};

  /* Get the top level image dimensions */
  const int top_num_pixels[2] = {
      ceil(c->top->image_data.padded_width[0] / image_data->pixel_size[0]),
      ceil(c->top->image_data.padded_width[1] / image_data->pixel_size[1])};

  /* Get the image dimensions for this sub cell */
  const int num_pixels[2] = {
      ceil(c->image_data.padded_width[0] / image_data->pixel_size[0]),
      ceil(c->image_data.padded_width[1] / image_data->pixel_size[1])};

  /* Get the top level location */
  const double xloc = c->top->image_data.padded_loc[0];
  const double yloc = c->top->image_data.padded_loc[1];

  /* Get the pixel location of this sub cell in the top level cells image. */
  const int pid = (int)((c->image_data.padded_loc[0] - xloc) / pixel_size[0]);
  const int pjd = (int)((c->image_data.padded_loc[1] - yloc) / pixel_size[1]);

  /* Get the sub cell image to add from. */
  double *sub_image = c->image_data.images[iimage];

  /* Loop over the pixels in the image. */
  for (int j = 0; j < num_pixels[0]; j++) {
    for (int k = 0; k < num_pixels[1]; k++) {
      /* Get the pixel location. */
      int x = pid + j;
      int y = pjd + k;

      /* Check if this pixel is in the image. */
      if (x < 0 || x >= top_num_pixels[0] || y < 0 || y >= top_num_pixels[1]) {
#ifdef SWIFT_DEBUG_CHECKS
        error("Pixel out of bounds: %d %d %d %d", x, y, top_num_pixels[0],
              top_num_pixels[1]);
#endif
        continue;
      }

      /* Add the pixel value to the top level image. */
      top_image[y + x * top_num_pixels[1]] += sub_image[k + j * num_pixels[1]];
    }
  }
}

void runner_do_imaging_collect(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  /* At the top level we need the images to collect into, this only the first
   * time round but we just do it at the point we need it.  */
  if (c->image_data.images == NULL && cell_init_images(c, r->e) != 0) {
    error("Failed to allocate memory for the top level image data.");
    return;
  }

  /* Get the image data for convenience. */
  struct image_common_data *image_data = r->e->image_data;
  int nimages = image_data->num_images;
  int num_pixels[2] = {c->top->image_data.num_pixels[0],
                       c->top->image_data.num_pixels[1]};

  /* Loop over all images. */
  for (int i = 0; i < nimages; i++) {
    /* Allocate an image to collect into. */
    double *top_image_buff;
    if (swift_memalign(
            "top_image_buff", (void **)&top_image_buff, SWIFT_STRUCT_ALIGNMENT,
            c->top->image_data.num_pixels[0] *
                c->top->image_data.num_pixels[1] * sizeof(double)) != 0) {
      error("Failed to allocate memory for the top level image.");
      return;
    }

    /* Zero the cell image buffer. */
    bzero(top_image_buff, c->top->image_data.num_pixels[0] *
                              c->top->image_data.num_pixels[1] *
                              sizeof(double));

    /* Collect the images from the sub cells. */
    runner_do_recursive_imaging_collect(r, c, top_image_buff, i);

    /* Copy over this image to the top level image. */
    double *top_image = c->top->image_data.images[i];
    for (int j = 0; j < num_pixels[0]; j++) {
      for (int k = 0; k < num_pixels[1]; k++) {
        int ipix = k + j * num_pixels[1];
        top_image[ipix] = top_image_buff[ipix];
      }
    }

    /* Free the top level image buffer. */
    free(top_image_buff);
  }

  if (timer) {
    TIMER_TOC(timer_dograv_down);
  }
}
