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
#include "runner.h"
#include "timers.h"

/**
 * @brief Make a mass image for dark matter.
 *
 * We treat dark matter as a top hat kernel with a radius of the softening
 * length.
 *
 * @param image_data The image data structure.
 * @param iimage The index of the image.
 * @param ci The cell data structure.
 */
static void runner_do_gpart_mass_image(struct image_data *image_data,
                                       int iimage, struct cell *c) {
  /* Extract the gparts from the cell. */
  struct gpart *gparts = c->grav.parts;
  int gcount = c->grav.count;

  /* Get the padded location we'll use as a zero point. */
  double xloc = c->image_data.padded_loc[0];
  double yloc = c->image_data.padded_loc[1];

  /* Get the pixel size. */
  double pixel_size[2] = {image_data->pixel_size[0], image_data->pixel_size[1]};

  /* Compute the number of pixels along each axis. */
  const int num_pixels[2] = {
      ceil(c->image_data.padded_width[0] / pixel_size[0]),
      ceil(c->image_data.padded_width[1] / pixel_size[1])};

  /* Extract the cells image buffer. */
  double *image = c->image_data.images[iimage];

  /* Loop over all gparts. */
  for (int i = 0; i < gcount; i++) {
    struct gpart *gp = &gparts[i];

    /* What pixel is this gpart in? */
    int pid = (int)((gp->x[0] - xloc) / image_data->pixel_size[0]);
    int pjd = (int)((gp->x[1] - yloc) / image_data->pixel_size[1]);

    /* For dark matter we smooth based on the softening length. How many pixels
     * does this cover? */
    int nsoft = (int)ceil((gp->epsilon / image_data->pixel_size[0]));

    /* If nsoft is 0/1 then we don't need to smooth. */
    // if (nsoft <= 1) {
    image[pjd + pid * num_pixels[1]] += gp->mass;
    continue;
    // }

    /* The square of the softening length for convenience. */
    double soft2 = gp->epsilon * gp->epsilon;

    /* Calculate the normalisation for a spherical top hat kernel with a
     * radius of nsoft pixels. */
    double norm = 1.0 / (M_PI * soft2 * soft2);

    /* Loop over the pixels in the smoothing kernel. */
    for (int j = -nsoft; j <= nsoft; j++) {
      for (int k = -nsoft; k <= nsoft; k++) {
        /* Get the pixel location. */
        int x = pid + j;
        int y = pjd + k;

        /* Check if this pixel is in the image. */
        if (x < 0 || x >= num_pixels[0] || y < 0 || y >= num_pixels[1]) {
#ifdef SWIFT_DEBUG_CHECKS
          error("Pixel out of bounds: %d %d %d %d", x, y, num_pixels[0],
                num_pixels[1]);
#endif
          continue;
        }

        /* Calculate the distance from the gpart to this pixel. */
        double dx = (x * image_data->pixel_size[0]) - gp->x[0];
        double dy = (y * image_data->pixel_size[1]) - gp->x[1];
        double r2 = dx * dx + dy * dy;

        /* If this is within the smoothing radius, add to the pixel value. */
        if (r2 < soft2) {
          image[y + x * num_pixels[1]] += gp->mass * norm;
        }
      }
    }
  }

  /* Flag that we generated an image. */
  c->image_data.image_generated = 1;
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

    /* Choose the right function for the job. */
    switch (ptype) {
      case 0: /* Gas */
        error("Gas imaging not implemented yet.");
        // if (strcmp(field_name, "density") == 0) {
        //   runner_do_part_density_image(&image_data->images[i], ci);
        // } else if (strcmp(field_name, "mass") == 0) {
        //   runner_do_part_mass_image(&image_data->images[i], ci);
        // } else if (strcmp(field_name, "temperature") == 0) {
        //   runner_do_part_temperature_image(&image_data->images[i], ci);
        // } else {
        //   error("Unknown field name for gas: %s", field_name);
        // }
        break;
      case 1: /* Dark matter */
        if (strcmp(field_name, "mass") == 0) {
          runner_do_gpart_mass_image(&image_data->images[i], i, c);
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
        //   runner_do_spart_mass_image(&image_data->images[i], ci);
        // } else {
        //   error("Unknown field name for stars: %s", field_name);
        // }
        break;
      case 5: /* Black holes */
        error("Black holes imaging not implemented yet.");
        // if (strcmp(field_name, "mass") == 0) {
        //   runner_do_bpart_mass_image(&image_data->images[i], ci);
        // } else {
        //   error("Unknown field name for black holes: %s", field_name);
        // }
        break;
      default:
        error("Unknown particle type: %d", ptype);
    }
  }

  if (timer) TIMER_TOC(timer_dograv_down);
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

  /* Reset the image generated flag */
  c->image_data.image_generated = 0;
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
            "top_image", (void **)&top_image_buff, SWIFT_STRUCT_ALIGNMENT,
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
