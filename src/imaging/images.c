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

/* This files header. */
#include "images.h"

/* Standard includes. */
#include <stddef.h>
#include <stdlib.h>

/* STB image writer for writing PNGs. */
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/* Local includes. */
#include "colormaps.h"
#include "parser.h"
#include "tools.h"

/**
 * @brief Intialise the image data structure.
 *
 * @param parameter_file The parsed parameter file.
 * @param e The engine data structure.
 *
 */
void imaging_init(struct swift_params *parameter_file, struct engine *e) {

  /* First flag that we are doing imaging. */
  e->with_imaging = 1;

  /* Allocate the image data structure. */
  if (swift_memalign("images_common", (void **)&e->image_data,
                     SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct image_common_data)) != 0) {
    error("Failed to allocate memory for the image data structure.");
    return;
  }

  /* Get the image data for convenience. */
  struct image_common_data *image_data = e->image_data;

  /* Read the common data from the parameter file. */
  image_data->num_images =
      parser_get_param_int(parameter_file, "ImagesCommon:nimages");
  if (image_data->num_images <= 0) {
    error("Number of images must be greater than 0.");
    return;
  }
  image_data->xres = parser_get_opt_param_int(
      parameter_file, "ImagesCommon:x_resolution", 1080);
  image_data->yres = parser_get_opt_param_int(
      parameter_file, "ImagesCommon:y_resolution", 1080);
  image_data->slice =
      parser_get_opt_param_int(parameter_file, "ImagesCommon:slice", 0);
  image_data->slice_thickness = parser_get_opt_param_double(
      parameter_file, "ImagesCommon:slice_thickness", 5.0);
  parser_get_opt_param_string(parameter_file, "ImagesCommon:subdir",
                              image_data->output_dir, "images");

  /* Compute the pixel size, this is boxsize / resolution. */
  image_data->pixel_size[0] = e->s->dim[0] / image_data->xres;
  image_data->pixel_size[1] = e->s->dim[1] / image_data->yres;

  /* Does the output directory exist? If not, create it. */
  if (e->nodeID == 0) {
    safe_checkdir(image_data->output_dir, /*create=*/1);
  }

  /* Allocate the images data structs. */
  if (swift_memalign("images", (void **)&image_data->images,
                     SWIFT_STRUCT_ALIGNMENT,
                     image_data->num_images * sizeof(struct image_data)) != 0) {
    error("Failed to allocate memory for the images.");
    return;
  }

  /* Read the image data from the parameter file from each ImageX block. */
  for (int i = 0; i < image_data->num_images; i++) {
    char block_name[256];
    snprintf(block_name, sizeof(block_name), "Image%d", i);
    struct image_data *image = &image_data->images[i];

    /* Read the image data from the parameter file. */
    char param_name[256];
    snprintf(param_name, sizeof(param_name), "%s:basename", block_name);
    parser_get_param_string(parameter_file, param_name, image->base_name);

    snprintf(param_name, sizeof(param_name), "%s:field_name", block_name);
    parser_get_param_string(parameter_file, param_name, image->field_name);

    snprintf(param_name, sizeof(param_name), "%s:particle_type", block_name);
    image->particle_type = parser_get_param_int(parameter_file, param_name);

    snprintf(param_name, sizeof(param_name), "%s:subdir", block_name);
    parser_get_opt_param_string(parameter_file, param_name, image->output_dir,
                                image_data->output_dir);

    snprintf(param_name, sizeof(param_name), "%s:slice", block_name);
    image->slice =
        parser_get_opt_param_int(parameter_file, param_name, image_data->slice);

    snprintf(param_name, sizeof(param_name), "%s:slice_thickness", block_name);
    image->slice_thickness = parser_get_opt_param_double(
        parameter_file, param_name, image_data->slice_thickness);

    /* Initialise the frame counter. */
    image->frame_number = 0;

    /* Copy over the iamge diemensions. */
    image->xres = image_data->xres;
    image->yres = image_data->yres;
    image->pixel_size[0] = image_data->pixel_size[0];
    image->pixel_size[1] = image_data->pixel_size[1];
  }

  /* Report some information for the hell of it. */
  if (e->verbose) {
    message("Number of images: %d", image_data->num_images);
    message("Output directory: %s", image_data->output_dir);
    if (image_data->slice) {
      message("  Slice mode: %s", image_data->slice ? "yes" : "no");
    }
    for (int i = 0; i < image_data->num_images; i++) {
      struct image_data *image = &image_data->images[i];
      message("Image %d: %s", i, image->base_name);
      message("Particle type: %d", image->particle_type);
      message("  Field name: %s", image->field_name);
      message("  Resolution: %dx%d", image->xres, image->yres);
      message("  Pixel size: %f x %f [Internal Units]", image->pixel_size[0],
              image->pixel_size[1]);
    }
  }
}

/**
 * @brief Write a PNG whose pixels are mapped through an arbitrary colormap.
 *
 * @param filename    Output filename (e.g. "out.png")
 * @param data        Input data array of length width*height (double)
 * @param width       Image width
 * @param height      Image height
 * @param cmap        Colormap: an array of [N][3] uint8_t entries {R,G,B}
 * @param cmap_size   Number of entries in cmap (e.g. 256)
 */
static void imaging_write_colormap_png_min_max(const char *filename,
                                               const double *data, int width,
                                               int height,
                                               const uint8_t cmap[][3],
                                               size_t cmap_size) {
  // 1) Find data min/max
  double mn = data[0], mx = data[0];
  size_t npix = (size_t)width * height;
  for (size_t i = 1; i < npix; i++) {
    if (data[i] < mn)
      mn = data[i];
    if (data[i] > mx)
      mx = data[i];
  }
  double inv_range = (mx > mn) ? 1.0 / (mx - mn) : 0.0;

  // 2) Allocate an RGB buffer
  unsigned char *rgb = malloc(3 * npix);
  if (!rgb)
    return;

  // 3) Map each sample into [0..cmap_size-1] and look up RGB
  for (size_t i = 0; i < npix; i++) {
    double norm = (data[i] - mn) * inv_range;
    if (norm < 0.0)
      norm = 0.0;
    if (norm > 1.0)
      norm = 1.0;
    size_t idx = (size_t)(norm * (cmap_size - 1) + 0.5);
    rgb[3 * i + 0] = cmap[idx][0];
    rgb[3 * i + 1] = cmap[idx][1];
    rgb[3 * i + 2] = cmap[idx][2];
  }

  // 4) Write to PNG (3 channels, stride = 3*width)
  stbi_write_png(filename, width, height, 3, rgb, 3 * width);

  free(rgb);
}

/**
 * @brief Write a PNG whose pixels are mapped through an arbitrary colormap
 *        with an automatic ±3σ contrast stretch and optional gamma.
 *
 * @param filename      Output filename (e.g. "out.png")
 * @param data          Input data array of length width*height (double)
 * @param width         Image width
 * @param height        Image height
 * @param cmap          Colormap: an array of [N][3] uint8_t entries {R,G,B}
 * @param cmap_size     Number of entries in cmap (e.g. 256)
 * @param gamma         Gamma to apply after stretch (e.g. 1.0 = linear,
 *                      1.8–2.2 to brighten midtones)
 */
static void imaging_write_colormap_png_zscale(const char *filename,
                                              const double *data, int width,
                                              int height,
                                              const uint8_t cmap[][3],
                                              size_t cmap_size, double gamma) {
  size_t npix = (size_t)width * height;

  // 1) Compute mean & variance
  double sum = 0.0, sum2 = 0.0;
  for (size_t i = 0; i < npix; ++i) {
    double v = data[i];
    sum += v;
    sum2 += v * v;
  }
  double mean = sum / npix;
  double var = sum2 / npix - mean * mean;
  double sigma = (var > 0.0 ? sqrt(var) : 0.0);

  // 2) Define our stretch window at ±3σ
  double lo = mean - 3.0 * sigma;
  double hi = mean + 3.0 * sigma;

  // 3) Fallback for flat images
  if (hi <= lo) {
    lo = mean - 1e-3;
    hi = mean + 1e-3;
  }
  double inv_range = 1.0 / (hi - lo);

  // 4) Allocate RGB buffer
  unsigned char *rgb = malloc(3 * npix);
  if (!rgb)
    return;

  // 5) Map each sample → [0..1], apply γ, then colormap index
  for (size_t i = 0; i < npix; ++i) {
    // linear stretch to [0,1]
    double norm = (data[i] - lo) * inv_range;
    if (norm < 0.0)
      norm = 0.0;
    else if (norm > 1.0)
      norm = 1.0;

    // optional gamma
    if (gamma != 1.0) {
      norm = pow(norm, 1.0 / gamma);
    }

    // lookup index
    size_t idx = (size_t)(norm * (cmap_size - 1) + 0.5);

    rgb[3 * i + 0] = cmap[idx][0];
    rgb[3 * i + 1] = cmap[idx][1];
    rgb[3 * i + 2] = cmap[idx][2];
  }

  // 6) Write PNG (row-major, stride = 3*width)
  stbi_write_png(filename, width, height, 3, rgb, 3 * width);

  free(rgb);
}

static void imaging_combine_cell_images(struct space *s,
                                        struct image_common_data *image_data,
                                        int image_number, double *image_buff) {

  /* Get the cells ready to loop over them. */
  struct cell *cells = s->cells_top;
  int ncells = s->nr_cells;

  /* Loop over all cells. */
  for (int i = 0; i < ncells; i++) {
    /* Get the cell. */
    struct cell *c = &cells[i];

    /* Skip empty cells. */
    if (c->hydro.count == 0 || c->grav.count == 0) {
      continue;
    }

    /* Extract the image data from this cell. */
    double *cell_image = c->image_data.images[image_number];
    double padded_loc[2] = {c->image_data.padded_loc[0],
                            c->image_data.padded_loc[1]};
    int num_pixels[2] = {c->image_data.num_pixels[0],
                         c->image_data.num_pixels[1]};

    /* Where is this cell in the whole image? */
    int pid = (int)(padded_loc[0] / image_data->pixel_size[0]);
    int pjd = (int)(padded_loc[1] / image_data->pixel_size[1]);

    /* Loop over the pixels in the image. */
    for (int j = 0; j < num_pixels[0]; j++) {
      for (int k = 0; k < num_pixels[1]; k++) {
        /* Get the pixel location in the image. */
        int xloc = pid + j;
        int yloc = pjd + k;

        /* Get the pixel index. */
        int idx = yloc + xloc * image_data->yres;

        /* Check if this pixel is in the image. */
        if (xloc < 0 || xloc >= image_data->xres || yloc < 0 ||
            yloc >= image_data->yres) {
#ifdef SWIFT_DEBUG_CHECKS
          error("Pixel out of bounds: %d %d %d %d", xloc, yloc,
                image_data->xres, image_data->yres);
#endif
          continue;
        }

        /* Check if this pixel is in the image. */
        if (k < 0 || k >= num_pixels[0] || j < 0 || j >= num_pixels[1]) {
#ifdef SWIFT_DEBUG_CHECKS
          error("Pixel out of bounds: %d %d %d %d", xloc, yloc, num_pixels[0],
                num_pixels[1]);
#endif
          continue;
        }

        /* Add the cell image to the image buffer. */
        image_buff[idx] += cell_image[k + j * num_pixels[1]];
      }
    }
  }
}

static void imaging_write_image(struct space *s,
                                struct image_common_data *image_data,
                                int image_number) {
  /* Get the image data for convenience. */
  struct image_data *image = &image_data->images[image_number];

  /* Allocate an image buffer to collect each cells image into. */
  double *image_buff;
  if (swift_memalign("image_buff", (void **)&image_buff, SWIFT_STRUCT_ALIGNMENT,
                     image_data->xres * image_data->yres * sizeof(double)) !=
      0) {
    error("Failed to allocate memory for the image buffer.");
    return;
  }
  bzero(image_buff, image_data->xres * image_data->yres * sizeof(double));

  /* Combine the cell images ready to be written out. */
  imaging_combine_cell_images(s, image_data, image_number, image_buff);

  /* Create the filename. */
  char filename[256];
  snprintf(filename, sizeof(filename), "%s/%s_%d.png", image_data->output_dir,
           image->base_name, image->frame_number);

  /* Write the image as an RGB PNG. */
  imaging_write_colormap_png_min_max(filename, image_buff, image_data->xres,
                                     image_data->yres, plasma_colormap,
                                     plasma_colormap_size);

  /* Free the image buffer. */
  free(image_buff);

  if (s->e->verbose) {
    message("Wrote image to %s", filename);
  }

  /* Increment the frame number. */
  image->frame_number++;
}

void imaging_write_images(struct engine *e) {
  /* Get the image data for convenience. */
  struct image_common_data *image_data = e->image_data;

  /* Loop over all images. */
  for (int i = 0; i < image_data->num_images; i++) {
    imaging_write_image(e->s, image_data, i);
  }
}
