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

  /* Get the image resolution. */
  image_data->xres =
      parser_get_opt_param_int(parameter_file, "Imaging:x_resolution", 1080);
  image_data->yres =
      parser_get_opt_param_int(parameter_file, "Imaging:y_resolution", 1080);

  /* Initialise the frame number at 0. */
  image_data->frame_number = 0;

  /* Where are we writing the images? */
  parser_get_opt_param_string(parameter_file, "Imaging:subdir",
                              image_data->output_dir, "images");

  /* Basename for the images. */
  parser_get_opt_param_string(parameter_file, "Imaging:basename",
                              image_data->base_name, "image");

  /* Does the output directory exist? If not, create it. */
  if (nodeID == 0) {
    safe_checkdir(image_data->output_dir, /*create=*/1);
  }

  /* Angular field of view in radians. */
  image_data->fov_angle[0] = M_PI / 3.0;  // 60 degrees
  image_data->fov_angle[1] = M_PI / 3.0;  // 60 degrees

  /* Camera distance from centre of the image. */
  image_data->camera_distance =
      parser_get_param_double(parameter_file, "Imaging:distance");

  /* Intialise the projected kernel table. */
  image_data->projected_kernel_table = (struct projected_kernel_table *)malloc(
      sizeof(struct projected_kernel_table));
  projected_kernel_init(image_data->projected_kernel_table);

  /* Report some information for the hell of it. */
  if (verbose) {
    message("Output directory: %s", image_data->output_dir);
    message("Base name: %s", image_data->base_name);
    message("Image resolution: %dx%d", image_data->xres, image_data->yres);
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

/**
 * @brief A mapper function to create all the FOS images at once for each cell.
 *
 * @param image_data The image data structure.
 * @param c The cell to process.
 */
void imaging_cell_mapper(void *map_data, int num_elements, void *extra_data) {

  /* Unpack the data we have been given. */
  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Get the iamge data we will need. */
  struct image_common_data *image_data = (struct image_common_data *)map_data;
  double image_centre[3] = {dim[0] * 0.5, dim[1] * 0.5, dim[2] * 0.5};
  double camera_distance = image_data->camera_distance;
  int xres = image_data->xres;
  int yres = image_data->yres;
  int angular_fov[2] = {image_data->fov_angle[0], image_data->fov_angle[1]};

  /* Scale factor: image plane spans [-scale_x, scale_x] and [-scale_y,
   * scale_y] */
  double scale_x = camera_distance * tan(angular_fov[0] / 2.0);
  double scale_y = camera_distance * tan(angular_fov[1] / 2.0);

  /* Threadpool id of current thread. */
  short int tpid = threadpool_gettid();

  /* Get each of the images for this thread. */
  double *dm_image = image_common_data->dm_images[tpid];
  double *gas_image = image_common_data->gas_images[tpid];
  double *stars_image = image_common_data->stars_images[tpid];
  double *gas_temp_image = image_common_data->gas_temp_images[tpid];

  /* Loop over the cells we have been given. */
  for (int i = 0; i < num_elements; i++) {

    /* Get the cell index. */
    int cid = (size_t)(map_data) + i;

    /* Get the cell. */
    struct cell *c = &cells[cid];

    /* Skip empty cells. */
    if (c->hydro.count == 0 && c->grav.count == 0) {
      continue;
    }

    /* Loop over all dark matter particles in this cell. */
    for (int j = 0; j < c->grav.count; j++) {
      /* Get the dark matter particle. */
      struct gpart *gp = &c->grav.parts[j];

      /* Skip if not a dark matter particle. */
      if (gp->type != swift_type_dark_matter) {
        continue;
      }

      /* Get the relative position of the particle to the image centre. */
      double pos[3] = {gp->pos[0] - image_centre[0],
                       gp->pos[1] - image_centre[1],
                       gp->pos[2] - image_centre[2]};

      /* Project the particle onto the image plane in angular coordinates at
       * the camera distance. */
      double theta = atan2(pos[1], pos[0]);
      double phi = atan2(pos[2], sqrt(pos[0] * pos[0] + pos[1] * pos[1]));
      double r = camera_distance * tan(phi);
      double x = r * cos(theta);
      double y = r * sin(theta);

      /* Normalized coordinates in [-1, 1] */
      double nx = x / scale_x;
      double ny = y / scale_y;

      /* Convert to pixel coordinates */
      int ix = (int)((nx + 1.0) * 0.5 * xres);
      int iy = (int)((ny + 1.0) * 0.5 * yres);

      /* Bounds check and increment image */
      if (ix >= 0 && ix < xres && iy >= 0 && iy < yres) {
        int index = iy * xres + ix;
        dm_image[index] += gp->mass;
      }
    }

    /* Loop over all gas particles in this cell. */
    for (int j = 0; j < c->hydro.count; j++) {
      /* Get the gas particle. */
      struct part *p = &c->hydro.parts[j];

      /* Get the relative position of the particle to the image centre. */
      double pos[3] = {p->pos[0] - image_centre[0], p->pos[1] - image_centre[1],
                       p->pos[2] - image_centre[2]};

      /* Project the particle onto the image plane in angular coordinates at
       * the camera distance. */
      double theta = atan2(pos[1], pos[0]);
      double phi = atan2(pos[2], sqrt(pos[0] * pos[0] + pos[1] * pos[1]));
      double r = camera_distance * tan(phi);
      double x = r * cos(theta);
      double y = r * sin(theta);

      /* Normalized coordinates in [-1, 1] */
      double nx = x / scale_x;
      double ny = y / scale_y;

      /* Convert to pixel coordinates */
      int ix = (int)((nx + 1.0) * 0.5 * xres);
      int iy = (int)((ny + 1.0) * 0.5 * yres);

      /* Bounds check and increment image */
      if (ix >= 0 && ix < xres && iy >= 0 && iy < yres) {
        int index = iy * xres + ix;
        gas_image[index] += p->mass;
        gas_temp_image[index] += p->cooling_data.subgrid_temp * p->mass;
      }
    }

    /* Loop over all star particles in this cell. */
    for (int j = 0; j < c->stars.count; j++) {
      /* Get the star particle. */
      struct spart *sp = &c->stars.parts[j];

      /* Get the relative position of the particle to the image centre. */
      double pos[3] = {sp->pos[0] - image_centre[0],
                       sp->pos[1] - image_centre[1],
                       sp->pos[2] - image_centre[2]};

      /* Project the particle onto the image plane in angular coordinates at
       * the camera distance. */
      double theta = atan2(pos[1], pos[0]);
      double phi = atan2(pos[2], sqrt(pos[0] * pos[0] + pos[1] * pos[1]));
      double r = camera_distance * tan(phi);
      double x = r * cos(theta);
      double y = r * sin(theta);

      /* Normalized coordinates in [-1, 1] */
      double nx = x / scale_x;
      double ny = y / scale_y;

      /* Convert to pixel coordinates */
      int ix = (int)((nx + 1.0) * 0.5 * xres);
      int iy = (int)((ny + 1.0) * 0.5 * yres);

      /* Bounds check and increment image */
      if (ix >= 0 && ix < xres && iy >= 0 && iy < yres) {
        int index = iy * xres + ix;
        stars_image[index] += sp->mass;
      }
    }
  }
}

void imaging_allocate_threadimages_mapper(void *map_data, int num_elements,
                                          void *extra_data) {
  /* Get the threadpool id of the current thread. */
  short int tpid = threadpool_gettid();

  /* Get this threads image data. */
  struct image_common_data *image_data = (struct image_common_data *)extra_data;

  /* Allocate the images for this thread. */
  if (swift_memalign("dm_images", (void **)&image_data->dm_images[tpid],
                     SWIFT_STRUCT_ALIGNMENT,
                     image_data->xres * image_data->yres * sizeof(double)) !=
      0) {
    error("Failed to allocate memory for the dark matter images.");
    return;
  }
  if (swift_memalign("gas_images", (void **)&image_data->gas_images[tpid],
                     SWIFT_STRUCT_ALIGNMENT,
                     image_data->xres * image_data->yres * sizeof(double)) !=
      0) {
    error("Failed to allocate memory for the gas images.");
    return;
  }
  if (swift_memalign("stars_images", (void **)&image_data->stars_images[tpid],
                     SWIFT_STRUCT_ALIGNMENT,
                     image_data->xres * image_data->yres * sizeof(double)) !=
      0) {
    error("Failed to allocate memory for the star images.");
    return;
  }
  if (swift_memalign(
          "gas_temp_images", (void **)&image_data->gas_temp_images[tpid],
          SWIFT_STRUCT_ALIGNMENT,
          image_data->xres * image_data->yres * sizeof(double)) != 0) {
    error("Failed to allocate memory for the gas temperature images.");
    return;
  }

  /* Zero the images. */
  bzero(image_data->dm_images[tpid],
        image_data->xres * image_data->yres * sizeof(double));
  bzero(image_data->gas_images[tpid],
        image_data->xres * image_data->yres * sizeof(double));
  bzero(image_data->stars_images[tpid],
        image_data->xres * image_data->yres * sizeof(double));
  bzero(image_data->gas_temp_images[tpid],
        image_data->xres * image_data->yres * sizeof(double));
}

/**
 * @brief Write images to HDF5 format.
 *
 * @param e The engine data structure.
 */
void imaging_write_images_hdf5(double *dm_image, double *gas_image,
                               double *stars_image, double *gas_temp_image,
                               int xres, int yres, const char *output_dir,
                               const char *base_name, int frame_number) {
  /* Create the output filename. */
  char filename[256];
  snprintf(filename, sizeof(filename), "%s/%s_%d.hdf5", output_dir,
           frame_number);

  /* Open the HDF5 file for writing. */
  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) {
    error("Failed to create HDF5 file %s.", filename);
    return;
  }

  /* Create the datasets for each image. */
  hsize_t dims[2] = {xres, yres};
  hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

  hid_t dm_dataset_id =
      H5Dcreate(file_id, "dark_matter", H5T_NATIVE_DOUBLE, dataspace_id,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t gas_dataset_id =
      H5Dcreate(file_id, "gas", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
  hid_t stars_dataset_id =
      H5Dcreate(file_id, "stars", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
  hid_t gas_temp_dataset_id =
      H5Dcreate(file_id, "gas_temperature", H5T_NATIVE_DOUBLE, dataspace_id,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Write the data to the datasets. */
  if (H5Dwrite(dm_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace_id,
               H5P_DEFAULT, dm_image) < 0) {
    error("Failed to write dark matter image to %s.", filename);
  }
  if (H5Dwrite(gas_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace_id,
               H5P_DEFAULT, gas_image) < 0) {
    error("Failed to write gas image to %s.", filename);
  }
  if (H5Dwrite(stars_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace_id,
               H5P_DEFAULT, stars_image) < 0) {
    error("Failed to write stars image to %s.", filename);
  }
  if (H5Dwrite(gas_temp_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace_id,
               H5P_DEFAULT, gas_temp_image) < 0) {
    error("Failed to write gas temperature image to %s.", filename);
  }

  /* Close the datasets and dataspace. */
  H5Dclose(dm_dataset_id);
  H5Dclose(gas_dataset_id);
  H5Dclose(stars_dataset_id);
  H5Dclose(gas_temp_dataset_id);
  H5Sclose(dataspace_id);

  /* Close the file. */
  H5Fclose(file_id);
  if (e->verbose) {
    message("Wrote images to %s", filename);
  }
}

/**
 * @brief Run the imaging loop all at once creating all the images at once.
 *
 * @param s The space structure containing the engine and cells.
 */
void imaging_compute_angular_images(struct space *s) {

  ticks tic = getticks();

  /* Unpack things we will need. */
  struct engine *e = s->e;
  struct image_common_data *image_data = e->image_data;

  /* Allocate the image buffers for each thread. */
  if (swift_memalign("dm_images", (void **)&image_data->dm_images,
                     SWIFT_STRUCT_ALIGNMENT,
                     e->threadpool.nthreads * sizeof(double *)) != 0) {
    error("Failed to allocate memory for the dark matter images.");
    return;
  }
  if (swift_memalign("gas_images", (void **)&image_data->gas_images,
                     SWIFT_STRUCT_ALIGNMENT,
                     e->threadpool.nthreads * sizeof(double *)) != 0) {
    error("Failed to allocate memory for the gas images.");
    return;
  }
  if (swift_memalign("stars_images", (void **)&image_data->stars_images,
                     SWIFT_STRUCT_ALIGNMENT,
                     e->threadpool.nthreads * sizeof(double *)) != 0) {
    error("Failed to allocate memory for the star images.");
    return;
  }
  if (swift_memalign("gas_temp_images", (void **)&image_data->gas_temp_images,
                     SWIFT_STRUCT_ALIGNMENT,
                     e->threadpool.nthreads * sizeof(double *)) != 0) {
    error("Failed to allocate memory for the gas temperature images.");
    return;
  }

  /* Now allocate the actual image buffers for each thread. */
  threadpool_map(&e->threadpool, imaging_allocate_threadimages_mapper, NULL,
                 e->threadpool.nthreads, 1, 1, image_data);

  /* Now we can calculate the images. */
  threadpool_map(&e->threadpool, imaging_cell_mapper, NULL, s->nr_cells, 1,
                 threadpool_auto_chunk_size, e);

  /* We now have all the thread images filled with the data. We now need to
   * reduce them into singular images. */
  double *dm_image = image_data->dm_images[0];
  double *gas_image = image_data->gas_images[0];
  double *stars_image = image_data->stars_images[0];
  double *gas_temp_image = image_data->gas_temp_images[0];
  for (int i = 1; i < e->threadpool.nthreads; i++) {
    /* Add the images from each thread together. */
    for (int j = 0; j < image_data->xres * image_data->yres; j++) {
      dm_image[j] += image_data->dm_images[i][j];
      gas_image[j] += image_data->gas_images[i][j];
      stars_image[j] += image_data->stars_images[i][j];
      gas_temp_image[j] += image_data->gas_temp_images[i][j];
    }
  }

  /* Final thing to do is divide out the mass weighting from the gas
   * temperature image. */
  for (int i = 0; i < image_data->xres * image_data->yres; i++) {
    if (gas_image[i] > 0.0) {
      gas_temp_image[i] /= gas_image[i];
    } else {
      gas_temp_image[i] = 0.0;
    }
  }

  /* Now we can write the images out into HDF5 files. */
  imaging_write_images_hdf5(dm_image, gas_image, stars_image, gas_temp_image,
                            image_data->xres, image_data->yres,
                            image_data->output_dir, image_data->base_name,
                            image_data->frame_number);

  /* Increment the frame number for the next time we write images. */
  image_data->frame_number++;

  /* Free the thread images. */
  for (int i = 0; i < e->threadpool.nthreads; i++) {
    free(image_data->dm_images[i]);
    free(image_data->gas_images[i]);
    free(image_data->stars_images[i]);
    free(image_data->gas_temp_images[i]);
  }
  free(image_data->dm_images);
  free(image_data->gas_images);
  free(image_data->stars_images);
  free(image_data->gas_temp_images);
  image_data->dm_images = NULL;
  image_data->gas_images = NULL;
  image_data->stars_images = NULL;
  image_data->gas_temp_images = NULL;

  if (e->verbose) {
    message("Computed angular images in %.3f seconds.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
  }
}
