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
#include "lightcone/projected_kernel.h"
#include "parser.h"
#include "tools.h"

/**
 * @brief Intialise the image data structure.
 *
 * @param parameter_file The parsed parameter file.
 * @param e The engine data structure.
 *
 */
void imaging_init(struct image_common_data *image_data,
                  struct swift_params *parameter_file, const int verbose,
                  const double dim[3], const int nodeID) {

  /* Read the common data from the parameter file. */

  /* How many images are we going to create? */
  image_data->num_images =
      parser_get_param_int(parameter_file, "ImagesCommon:nimages");
  if (image_data->num_images <= 0) {
    error("Number of images must be greater than 0.");
    return;
  }

  /* Get the image resolution. */
  image_data->xres = parser_get_opt_param_int(
      parameter_file, "ImagesCommon:x_resolution", 1080);
  image_data->yres = parser_get_opt_param_int(
      parameter_file, "ImagesCommon:y_resolution", 1080);

  /* Are we only doing a slice? */
  image_data->slice =
      parser_get_opt_param_int(parameter_file, "ImagesCommon:slice", 0);
  image_data->slice_thickness = parser_get_opt_param_double(
      parameter_file, "ImagesCommon:slice_thickness", 5.0);

  /* Where are we writing the images? */
  parser_get_opt_param_string(parameter_file, "ImagesCommon:subdir",
                              image_data->output_dir, "images");

  /* PNGs or raw arrays? */
  image_data->write_pngs =
      parser_get_opt_param_int(parameter_file, "ImagesCommon:write_pngs", 1);
  image_data->write_raw_arrays = !image_data->write_pngs;

  /* Are we doing a subvolume? */
  image_data->subvolume =
      parser_get_opt_param_int(parameter_file, "ImagesCommon:subvolume", 0);
  image_data->subvolume_centre[0] = 0.0;
  image_data->subvolume_centre[1] = 0.0;
  image_data->subvolume_centre[2] = 0.0;
  parser_get_opt_param_double_array(parameter_file, "ImagesCommon:centre", 3,
                                    image_data->subvolume_centre);
  if (image_data->subvolume) {
    /* If we are doing a subvolume, get the size of the field of view. */
    parser_get_opt_param_double_array(parameter_file, "ImagesCommon:fov", 3,
                                      image_data->fov);
    if (image_data->fov[0] <= 0 || image_data->fov[1] <= 0 ||
        image_data->fov[2] <= 0) {
      error("Field of view must be positive in all dimensions.");
      return;
    }
  } else {
    /* If we are not doing a subvolume, set the FOV to the box size. */
    image_data->fov[0] = dim[0];
    image_data->fov[1] = dim[1];
    image_data->fov[2] = dim[2];
  }

  /* Set the origin of the image (this is either 0.0, 0.0, 0.0 or the edge of
   * the subvolume). */
  if (image_data->subvolume) {
    image_data->origin[0] =
        image_data->subvolume_centre[0] - image_data->fov[0] / 2.0;
    image_data->origin[1] =
        image_data->subvolume_centre[1] - image_data->fov[1] / 2.0;
    image_data->origin[2] =
        image_data->subvolume_centre[2] - image_data->fov[2] / 2.0;
  } else {
    image_data->origin[0] = 0.0;
    image_data->origin[1] = 0.0;
    image_data->origin[2] = 0.0;
  }

  /* Compute the pixel size, this is boxsize / resolution. */
  image_data->pixel_size[0] = image_data->fov[0] / image_data->xres;
  image_data->pixel_size[1] = image_data->fov[1] / image_data->yres;

  /* Does the output directory exist? If not, create it. */
  if (nodeID == 0) {
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

    /* Check that the output directory exists. */
    if (nodeID == 0) {
      safe_checkdir(image->output_dir, /*create=*/1);
    }

    snprintf(param_name, sizeof(param_name), "%s:slice", block_name);
    image->slice =
        parser_get_opt_param_int(parameter_file, param_name, image_data->slice);

    snprintf(param_name, sizeof(param_name), "%s:slice_thickness", block_name);
    image->slice_thickness = parser_get_opt_param_double(
        parameter_file, param_name, image_data->slice_thickness);

    /* Is this image weighted by another image? */
    snprintf(param_name, sizeof(param_name), "%s:weight_by", block_name);
    image->weight_by = parser_get_opt_param_int(parameter_file, param_name, -1);

    /* Initialise the frame counter. */
    image->frame_number = 0;

    /* Copy over the image diemensions. */
    image->xres = image_data->xres;
    image->yres = image_data->yres;
    image->pixel_size[0] = image_data->pixel_size[0];
    image->pixel_size[1] = image_data->pixel_size[1];
  }

  /* Ensure any weighting is being used on the same particle types and that
   * the image being weighted exists. */
  for (int i = 0; i < image_data->num_images; i++) {
    struct image_data *image = &image_data->images[i];
    if (image->weight_by >= 0) {
      if (image->weight_by >= image_data->num_images) {
        error("Image %d is weighted by image %d, which does not exist.",
              image->weight_by, i);
      }
      if (image->particle_type !=
          image_data->images[image->weight_by].particle_type) {
        error("Image %d is weighted by image %d, but they are different "
              "particle types (%d vs %d).",
              image->weight_by, i, image->particle_type,
              image_data->images[image->weight_by].particle_type);
      }
    }
  }

  /* Intialise the projected kernel table. */
  image_data->projected_kernel_table = (struct projected_kernel_table *)malloc(
      sizeof(struct projected_kernel_table));
  projected_kernel_init(image_data->projected_kernel_table);

  /* Report some information for the hell of it. */
  if (verbose) {
    message("Number of images: %d", image_data->num_images);
    message("Output directory: %s", image_data->output_dir);
    if (image_data->slice) {
      message("  Slice thickness: %g", image_data->slice_thickness);
    }
    message("Image resolution: %dx%d", image_data->xres, image_data->yres);
    message("Pixel size: %g x %g", image_data->pixel_size[0],
            image_data->pixel_size[1]);
    message("Writing PNGs: %s", image_data->write_pngs ? "yes" : "no");
    message("Writing raw arrays: %s",
            image_data->write_raw_arrays ? "yes" : "no");
    for (int i = 0; i < image_data->num_images; i++) {
      struct image_data *image = &image_data->images[i];
      message("Image %d: %s", i, image->base_name);
      message(" - Particle type: %d", image->particle_type);
      message(" - Field name: %s", image->field_name);
      if (image->output_dir != image_data->output_dir) {
        message(" - Output directory: %s", image->output_dir);
      }
      if (!image_data->slice && image->slice) {
        message(" - Slice thickness: %g", image->slice_thickness);
      }
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

    // /* Skip cells that are not in the image. */
    // if (!imaging_cell_overlaps_fov(image_data, c)) {
    //   continue;
    // }

    /* Extract the image data from this cell. */
    double *cell_image = c->image_data.images[image_number];
    double padded_loc[2] = {c->image_data.padded_loc[0] - image_data->origin[0],
                            c->image_data.padded_loc[1] -
                                image_data->origin[1]};
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

        /* Apply periodic boundary conditions if not doing a subvolume. */
        if (!image_data->subvolume) {
          xloc = (xloc + image_data->xres) % image_data->xres;
          yloc = (yloc + image_data->yres) % image_data->yres;
        } else {
          /* If we are doing a subvolume, we need to check if the pixel is
           * within the subvolume. */
          if (xloc < 0 || xloc >= image_data->xres || yloc < 0 ||
              yloc >= image_data->yres) {
            continue; // Skip pixels outside the subvolume.
          }
        }
        xloc = (xloc + image_data->xres) % image_data->xres;
        yloc = (yloc + image_data->yres) % image_data->yres;

        /* Get the pixel index. */
        int idx = yloc + xloc * image_data->yres;

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

static void imaging_write_image_raw(const char *filename,
                                    struct image_common_data *image_data,
                                    struct image_data *image,
                                    double *image_buff) {

  /* Write the image as a raw array. */
  FILE *fp = fopen(filename, "wb");
  if (fp == NULL) {
    error("Failed to open file %s for writing.", filename);
    return;
  }

  /* Write the image data. */
  size_t written = fwrite(image_buff, sizeof(double),
                          image_data->xres * image_data->yres, fp);
  if (written != (size_t)(image_data->xres * image_data->yres)) {
    error("Failed to write all data to file %s.", filename);
    fclose(fp);
    return;
  }

  /* Close the file. */
  fclose(fp);
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

  /* Are we writing a PNG or a raw array? */
  char filename[256];
  if (image_data->write_pngs) {
    /* Create the filename. */
    snprintf(filename, sizeof(filename), "%s/%s_%d.png", image->output_dir,
             image->base_name, image->frame_number);

    /* Write the image as an RGB PNG. */
    imaging_write_colormap_png_min_max(filename, image_buff, image_data->xres,
                                       image_data->yres, plasma_colormap,
                                       plasma_colormap_size);
  } else if (image_data->write_raw_arrays) {
    /* Create the filename. */
    snprintf(filename, sizeof(filename), "%s/%s_%d.dat", image->output_dir,
             image->base_name, image->frame_number);
    imaging_write_image_raw(filename, image_data, image, image_buff);
  }

  /* Free the image buffer. */
  free(image_buff);

  if (s->e->verbose) {
    message("Wrote image to %s", filename);
  }

  /* Increment the frame number. */
  image->frame_number++;
}

static void imaging_write_weighted_image(struct space *s,
                                         struct image_common_data *image_data,
                                         int image_number) {
  /* Get the image data for convenience. */
  struct image_data *image = &image_data->images[image_number];

  /* Get the image for weighting too. */
  struct image_data *weight_image = &image_data->images[image->weight_by];

  /* Allocate an image buffer to collect each cells image into. */
  double *image_buff;
  if (swift_memalign("image_buff", (void **)&image_buff, SWIFT_STRUCT_ALIGNMENT,
                     image_data->xres * image_data->yres * sizeof(double)) !=
      0) {
    error("Failed to allocate memory for the image buffer.");
    return;
  }
  bzero(image_buff, image_data->xres * image_data->yres * sizeof(double));

  /* Allocate an image buffer to for the weights. */
  double *weight_buff;
  if (swift_memalign(
          "weight_buff", (void **)&weight_buff, SWIFT_STRUCT_ALIGNMENT,
          image_data->xres * image_data->yres * sizeof(double)) != 0) {
    error("Failed to allocate memory for the weight buffer.");
    return;
  }
  bzero(weight_buff, image_data->xres * image_data->yres * sizeof(double));

  /* Combine the cell images ready to be written out. */
  imaging_combine_cell_images(s, image_data, image_number, image_buff);
  imaging_combine_cell_images(s, image_data, image->weight_by, weight_buff);

  /* Loop over the pixels dividing out the weights. */
  size_t npix = (size_t)image_data->xres * image_data->yres;
  for (size_t i = 0; i < npix; i++) {
    if (weight_buff[i] > 0.0) {
      image_buff[i] /= weight_buff[i];
    } else {
      image_buff[i] = 0.0;
    }
  }

  /* We're done with the weights now. */
  free(weight_buff);

  /* Are we writing a PNG or a raw array? */
  char filename[256];
  if (image_data->write_pngs) {
    /* Create the filename. */
    snprintf(filename, sizeof(filename), "%s/%s_%d.png", image->output_dir,
             image->base_name, image->frame_number);

    /* Write the image as an RGB PNG. */
    imaging_write_colormap_png_min_max(filename, image_buff, image_data->xres,
                                       image_data->yres, plasma_colormap,
                                       plasma_colormap_size);
  } else if (image_data->write_raw_arrays) {
    /* Create the filename. */
    snprintf(filename, sizeof(filename), "%s/%s_%d.dat", image->output_dir,
             image->base_name, image->frame_number);
    imaging_write_image_raw(filename, image_data, image, image_buff);
  }

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

    /* Are we weighting? */
    if (image_data->images[i].weight_by >= 0) {

      /* Write the image as a weighted image (this involves dividing out the
       * weights). */
      imaging_write_weighted_image(e->s, image_data, i);

    } else {

      /* Write the image as a normal image. */
      imaging_write_image(e->s, image_data, i);
    }
  }

  /* Ok, we are done. Reset the flag for imaging so it can be reevaluated
   * next timestep. */
  e->imaging_this_timestep = 0;
}

/**
 * @brief Clean up the imaging data.
 *
 * @param image_data The image data to clean up.
 */
void imaging_clean(struct image_common_data *image_data) {
  if (image_data == NULL) {
    return;
  }

  /* Free the images. */
  if (image_data->images != NULL) {
    free(image_data->images);
    image_data->images = NULL;
  }

  /* Free the projected kernel table. */
  if (image_data->projected_kernel_table != NULL) {
    projected_kernel_clean(image_data->projected_kernel_table);
    free(image_data->projected_kernel_table);
    image_data->projected_kernel_table = NULL;
  }

  /* Free the image data structure itself. */
  free(image_data);
  image_data = NULL;
}

/**
 * @brief Test if a cell is within the FOV.
 *
 * @param image_data The image data structure.
 * @param c The cell to test.
 *
 * @return 1 if the cell is within the FOV, 0 otherwise.
 */
int imaging_cell_overlaps_fov(const struct image_common_data *image_data,
                              const struct cell *c) {
  /* If we are not doing a subvolume, all cells are in the FOV. */
  if (!image_data->subvolume) {
    return 1;
  }

  /* Padded cell boundaries */
  const double cell_min[3] = {c->image_data.padded_loc[0],
                              c->image_data.padded_loc[1],
                              c->image_data.padded_loc[2]};
  const double cell_max[3] = {
      c->image_data.padded_loc[0] + c->image_data.padded_width[0],
      c->image_data.padded_loc[1] + c->image_data.padded_width[1],
      c->image_data.padded_loc[2] + c->image_data.padded_width[2]};

  /* FOV boundaries */
  const double fov_min[3] = {image_data->origin[0], image_data->origin[1],
                             image_data->origin[2]};
  const double fov_max[3] = {
      image_data->origin[0] + image_data->fov[0],
      image_data->origin[1] + image_data->fov[1],
      image_data->origin[2] + image_data->fov[2],
  };

  /* Calculate overlap in each dimension */
  const double overlap_x =
      fmin(cell_max[0], fov_max[0]) - fmax(cell_min[0], fov_min[0]);
  const double overlap_y =
      fmin(cell_max[1], fov_max[1]) - fmax(cell_min[1], fov_min[1]);
  const double overlap_z =
      fmin(cell_max[2], fov_max[2]) - fmax(cell_min[2], fov_min[2]);

  /* Check if overlap lengths are positive */
  return (overlap_x > 0.0 && overlap_y > 0.0 && overlap_z > 0.0);
}
