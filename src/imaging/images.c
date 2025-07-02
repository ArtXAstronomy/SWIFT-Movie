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

/* Local includes. */
#include "colormaps.h"
#include "lightcone/projected_kernel.h"
#include "parser.h"
#include "threadpool.h"
#include "tools.h"

/**
 * @brief Lazily allocate (if needed) and zero per-thread image buffers.
 *
 * This replaces both your old alloc‐mapper and zero‐mapper.
 */
void imaging_allocate_threadimages(struct engine *e) {

  /* Unpack useful data. */
  struct image_common_data *image_data = e->image_data;
  size_t npix = (size_t)image_data->xres * image_data->yres;
  size_t nbytes = npix * sizeof(double);

  /* Allocate each channel once */
  for (int tpid = 0; tpid < e->nr_threads; tpid++) {

    if (image_data->dm_images[tpid] == NULL) {
      if (swift_memalign("dm_images", (void **)&image_data->dm_images[tpid],
                         SWIFT_STRUCT_ALIGNMENT, nbytes) != 0) {
        error("Failed to alloc DM image for thread %d", tpid);
        return;
      }
    }
    if (image_data->gas_images[tpid] == NULL) {
      if (swift_memalign("gas_images", (void **)&image_data->gas_images[tpid],
                         SWIFT_STRUCT_ALIGNMENT, nbytes) != 0) {
        error("Failed to alloc gas image for thread %d", tpid);
        return;
      }
    }
    if (image_data->star_images[tpid] == NULL) {
      if (swift_memalign("star_images", (void **)&image_data->star_images[tpid],
                         SWIFT_STRUCT_ALIGNMENT, nbytes) != 0) {
        error("Failed to alloc star image for thread %d", tpid);
        return;
      }
    }
    if (image_data->gas_temp_images[tpid] == NULL) {
      if (swift_memalign("gas_temp_images",
                         (void **)&image_data->gas_temp_images[tpid],
                         SWIFT_STRUCT_ALIGNMENT, nbytes) != 0) {
        error("Failed to alloc gas temp image for thread %d", tpid);
        return;
      }
    }

    /* Zero the images. */
    bzero(image_data->dm_images[tpid], nbytes);
    bzero(image_data->gas_images[tpid], nbytes);
    bzero(image_data->star_images[tpid], nbytes);
    bzero(image_data->gas_temp_images[tpid], nbytes);
  }
}

/**
 * @brief Intialise the image data structure.
 *
 * @param parameter_file The parsed parameter file.
 * @param e The engine data structure.
 *
 */
void imaging_init(struct image_common_data *image_data,
                  struct swift_params *parameter_file, const int verbose,
                  const double dim[3], const int nodeID, const int nr_threads) {

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

  /* Create the staging directory for the images too if it does not exist. */
  char staging_dir[256];
  snprintf(staging_dir, sizeof(staging_dir), "%s_tmp", image_data->output_dir);
  if (nodeID == 0) {
    safe_checkdir(staging_dir, /*create=*/1);
  }

  /* Angular field of view in radians, here we need to account for whether the
   * image is square or not. */
  const double fixed_fov = M_PI / 3.0; /* 60° */
  double hfov, vfov;
  if (image_data->xres > image_data->yres) {
    hfov = fixed_fov;
    vfov = 2.0 * atan(tan(fixed_fov * 0.5) *
                      ((double)image_data->yres / image_data->xres));
  } else if (image_data->yres > image_data->xres) {
    vfov = fixed_fov;
    hfov = 2.0 * atan(tan(fixed_fov * 0.5) *
                      ((double)image_data->xres / image_data->yres));
  } else {
    /* square pixels: both axes = 60° */
    hfov = vfov = fixed_fov;
  }

  image_data->fov_angle[0] = hfov;
  image_data->fov_angle[1] = vfov;

  /* Get the camera position in spherical coordinates (R, theta, phi). */
  image_data->sphere_camera_position[0] =
      parser_get_opt_param_double(parameter_file, "Imaging:distance", 1.0);
  image_data->sphere_camera_position[1] =
      parser_get_opt_param_double(parameter_file, "Imaging:theta", M_PI / 4.0);
  image_data->sphere_camera_position[2] =
      parser_get_opt_param_double(parameter_file, "Imaging:phi", M_PI / 4.0);

  /* Intialise the projected kernel table. */
  image_data->projected_kernel_table = (struct projected_kernel_table *)malloc(
      sizeof(struct projected_kernel_table));
  projected_kernel_init(image_data->projected_kernel_table);

  /* Get the number of rotation frames. */
  image_data->nr_rotation_frames = parser_get_opt_param_int(
      parameter_file, "Imaging:nr_rotation_frames", 180);

  /* Allocate the image buffers for each thread. */
  if (swift_memalign("dm_images", (void **)&image_data->dm_images,
                     SWIFT_STRUCT_ALIGNMENT,
                     nr_threads * sizeof(double *)) != 0) {
    error("failed to allocate memory for the dark matter images.");
    return;
  }
  if (swift_memalign("gas_images", (void **)&image_data->gas_images,
                     SWIFT_STRUCT_ALIGNMENT,
                     nr_threads * sizeof(double *)) != 0) {
    error("failed to allocate memory for the gas images.");
    return;
  }
  if (swift_memalign("star_images", (void **)&image_data->star_images,
                     SWIFT_STRUCT_ALIGNMENT,
                     nr_threads * sizeof(double *)) != 0) {
    error("failed to allocate memory for the star images.");
    return;
  }
  if (swift_memalign("gas_temp_images", (void **)&image_data->gas_temp_images,
                     SWIFT_STRUCT_ALIGNMENT,
                     nr_threads * sizeof(double *)) != 0) {
    error("failed to allocate memory for the gas temperature images.");
    return;
  }

  /* Allocate the images for each thread. */
  for (int tpid = 0; tpid < nr_threads; tpid++) {
    image_data->dm_images[tpid] = NULL;
    image_data->gas_images[tpid] = NULL;
    image_data->star_images[tpid] = NULL;
    image_data->gas_temp_images[tpid] = NULL;
  }

  /* Report some information for the hell of it. */
  if (verbose) {
    message("Output directory: %s", image_data->output_dir);
    message("Base name: %s", image_data->base_name);
    message("Image resolution: %dx%d", image_data->xres, image_data->yres);
    message("Field of view: %.2f x %.2f radians", image_data->fov_angle[0],
            image_data->fov_angle[1]);
    message("Camera position (spherical): R=%.2f, theta=%.2f, phi=%.2f",
            image_data->sphere_camera_position[0],
            image_data->sphere_camera_position[1],
            image_data->sphere_camera_position[2]);
    message("Number of rotation frames: %d", image_data->nr_rotation_frames);
  }

  /* If gui_data.txt already exists in simulation directory, remove it. */
  if (nodeID == 0) {
    FILE *f = fopen("gui_data.txt", "w");
    fprintf(f, "  %6s %12s %12s %12s %12s %12s %12s %21s %12s\n", "step", "a",
            "z", "Nparts", "Ngparts", "Nsparts", "Nbparts", "Wallclock", "%");
    fclose(f);
  }
}

void imaging_write_images(struct engine *e) {

  /* Call the imaging function to generate and write the images. */
  struct space *s = e->s;
  imaging_compute_angular_images(s);

  /* Ok, we are done. Reset the flag for imaging so it can be reevaluated
   * next timestep. */
  e->imaging_this_timestep = 0;
}

/**
 * @brief Clean up the imaging data.
 *
 * @param image_data The image data to clean up.
 */
void imaging_clean(struct image_common_data *image_data, const int nr_threads) {
  if (image_data == NULL) {
    return;
  }

  /* Free the projected kernel table. */
  if (image_data->projected_kernel_table != NULL) {
    projected_kernel_clean(image_data->projected_kernel_table);
    free(image_data->projected_kernel_table);
    image_data->projected_kernel_table = NULL;
  }

  /* Free the thread images. */
  for (int tpid = 0; tpid < nr_threads; tpid++) {
    if (image_data->dm_images[tpid] != NULL) {
      free(image_data->dm_images[tpid]);
      image_data->dm_images[tpid] = NULL;
    }
    if (image_data->gas_images[tpid] != NULL) {
      free(image_data->gas_images[tpid]);
      image_data->gas_images[tpid] = NULL;
    }
    if (image_data->star_images[tpid] != NULL) {
      free(image_data->star_images[tpid]);
      image_data->star_images[tpid] = NULL;
    }
    if (image_data->gas_temp_images[tpid] != NULL) {
      free(image_data->gas_temp_images[tpid]);
      image_data->gas_temp_images[tpid] = NULL;
    }
  }
  free(image_data->dm_images);
  image_data->dm_images = NULL;
  free(image_data->gas_images);
  image_data->gas_images = NULL;
  free(image_data->star_images);
  image_data->star_images = NULL;
  free(image_data->gas_temp_images);

  /* Free the image data structure itself. */
  free(image_data);
  image_data = NULL;
}

/**
 * @brief A high-performance mapper to project & smooth particles into images
 *
 * Implements a single-pass SPH‐kernel deposition per particle, hoisting
 * invariants and minimizing work inside the tight inner loops.
 *
 * @param map_data     (size_t) start cell index for this thread
 * @param num_elements number of cells to process
 * @param extra_data   (struct engine*) containing image_data & space
 */
void imaging_cell_mapper(void *map_data, int num_elements, void *extra_data) {
  /* Thread, engine, and space setup */
  short tpid = threadpool_gettid();
  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  struct image_common_data *id = e->image_data;

  /* Image & FOV parameters */
  int xres = id->xres;
  int yres = id->yres;
  double fov_x = id->fov_angle[0];
  double fov_y = id->fov_angle[1];
  double R = id->sphere_camera_position[0];
  struct projected_kernel_table *kt = id->projected_kernel_table;
  double u_max = kt->u_max;
  double u_max_sq = u_max * u_max;

  /* Precompute angular→pixel scales */
  double max_x = tan(fov_x * 0.5);
  double max_y = tan(fov_y * 0.5);
  double px_per_rad_x = xres / (2.0 * max_x);
  double px_per_rad_y = yres / (2.0 * max_y);

  /* Build camera basis (forward, right, up) from spherical coords */
  double theta = id->sphere_camera_position[1];
  double phi = id->sphere_camera_position[2];
  double cam_rel[3] = {R * sin(theta) * cos(phi), R * sin(theta) * sin(phi),
                       R * cos(theta)};
  /* Forward = normalize(-cam_rel) */
  double F[3] = {-cam_rel[0], -cam_rel[1], -cam_rel[2]};
  double norm = sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
  F[0] /= norm;
  F[1] /= norm;
  F[2] /= norm;
  /* Adjust world‐up if nearly colinear */
  double world_up[3] = {0.0, 0.0, 1.0};
  if (fabs(F[0] * world_up[0] + F[1] * world_up[1] + F[2] * world_up[2]) >
      0.999) {
    world_up[0] = 0;
    world_up[1] = 1;
    world_up[2] = 0;
  }
  /* Right = normalize(cross(world_up, forward)) */
  double Rvec[3] = {world_up[1] * F[2] - world_up[2] * F[1],
                    world_up[2] * F[0] - world_up[0] * F[2],
                    world_up[0] * F[1] - world_up[1] * F[0]};
  norm = sqrt(Rvec[0] * Rvec[0] + Rvec[1] * Rvec[1] + Rvec[2] * Rvec[2]);
  Rvec[0] /= norm;
  Rvec[1] /= norm;
  Rvec[2] /= norm;
  /* Up = cross(forward, right) */
  double Uvec[3] = {F[1] * Rvec[2] - F[2] * Rvec[1],
                    F[2] * Rvec[0] - F[0] * Rvec[2],
                    F[0] * Rvec[1] - F[1] * Rvec[0]};

  /* Thread‐local image pointers */
  double *dm_img = id->dm_images[tpid];
  double *gas_img = id->gas_images[tpid];
  double *star_img = id->star_images[tpid];
  double *gtmp_img = id->gas_temp_images[tpid];

  /* Base cell index for this chunk */
  size_t base = (size_t)map_data;

  /* Main loop over cells assigned to this thread */
  for (int idx = 0; idx < num_elements; idx++) {
    int cid = (int)(base + idx);
    if (cid < 0 || cid >= s->nr_cells) continue;
    struct cell *c = &cells[cid];

    /* === Dark Matter === */
    for (int j = 0; j < c->grav.count; j++) {
      struct gpart *gp = &c->grav.parts[j];
      if (gp->type != swift_type_dark_matter) continue;

      /* World → camera → angular → pixel */
      double P[3] = {gp->x[0] - s->dim[0] * 0.5, gp->x[1] - s->dim[1] * 0.5,
                     gp->x[2] - s->dim[2] * 0.5};
      double V[3] = {P[0] - cam_rel[0], P[1] - cam_rel[1], P[2] - cam_rel[2]};
      double x_cam = V[0] * Rvec[0] + V[1] * Rvec[1] + V[2] * Rvec[2];
      double y_cam = V[0] * Uvec[0] + V[1] * Uvec[1] + V[2] * Uvec[2];
      double z_cam = V[0] * F[0] + V[1] * F[1] + V[2] * F[2];
      if (z_cam <= 1e-6) continue;
      double fx = (x_cam / z_cam / max_x * 0.5 + 0.5) * xres;
      double fy = (y_cam / z_cam / max_y * 0.5 + 0.5) * yres;

      /* Smoothing length in pixels & inv */
      double h_pix = (gp->epsilon * 4 / R) * px_per_rad_x;
      double inv_hpix = 1.0 / h_pix;

      /* Tiny‐kernel fallback */
      if (h_pix < 0.5) {
        int ix = (int)floor(fx + 0.5), iy = (int)floor(fy + 0.5);
        if (ix >= 0 && ix < xres && iy >= 0 && iy < yres)
          dm_img[iy * xres + ix] += gp->mass;
        continue;
      }

      /* Determine stencil radius */
      int fx_i = (int)floor(fx);
      int fy_i = (int)floor(fy);
      int delta = (int)ceil(u_max * h_pix);

      /* One‐pass deposition */
      for (int dj = -delta; dj <= delta; dj++) {
        int iy = fy_i + dj;
        if (iy < 0 || iy >= yres) continue;
        int row = iy * xres;
        double ry2 = pow(fabs(fy - (iy + 0.5)) * inv_hpix, 2);

        for (int di = -delta; di <= delta; di++) {
          int ix = fx_i + di;
          if (ix < 0 || ix >= xres) continue;
          double rx2 = pow(fabs(fx - (ix + 0.5)) * inv_hpix, 2);
          double u2 = rx2 + ry2;
          if (u2 >= u_max_sq) continue;
          double w = projected_kernel_eval(kt, sqrt(u2));
          dm_img[row + ix] += gp->mass * w;
        }
      }
    }

    /* === Gas === */
    for (int j = 0; j < c->hydro.count; j++) {
      struct part *p = &c->hydro.parts[j];

      double P[3] = {p->x[0] - s->dim[0] * 0.5, p->x[1] - s->dim[1] * 0.5,
                     p->x[2] - s->dim[2] * 0.5};
      double V[3] = {P[0] - cam_rel[0], P[1] - cam_rel[1], P[2] - cam_rel[2]};
      double x_cam = V[0] * Rvec[0] + V[1] * Rvec[1] + V[2] * Rvec[2];
      double y_cam = V[0] * Uvec[0] + V[1] * Uvec[1] + V[2] * Uvec[2];
      double z_cam = V[0] * F[0] + V[1] * F[1] + V[2] * F[2];
      if (z_cam <= 1e-6) continue;
      double fx = (x_cam / z_cam / max_x * 0.5 + 0.5) * xres;
      double fy = (y_cam / z_cam / max_y * 0.5 + 0.5) * yres;

      double h_pix = (p->h / R) * px_per_rad_x;
      double inv_hpix = 1.0 / h_pix;
      if (h_pix < 0.5) {
        int ix = (int)floor(fx + 0.5), iy = (int)floor(fy + 0.5);
        if (ix >= 0 && ix < xres && iy >= 0 && iy < yres) {
          int pix = iy * xres + ix;
          gas_img[pix] += p->mass;
          gtmp_img[pix] += p->mass * p->cooling_data.subgrid_temp;
        }
        continue;
      }

      int fx_i = (int)floor(fx);
      int fy_i = (int)floor(fy);
      int delta = (int)ceil(u_max * h_pix);

      for (int dj = -delta; dj <= delta; dj++) {
        int iy = fy_i + dj;
        if (iy < 0 || iy >= yres) continue;
        int row = iy * xres;
        double ry2 = pow(fabs(fy - (iy + 0.5)) * inv_hpix, 2);

        for (int di = -delta; di <= delta; di++) {
          int ix = fx_i + di;
          if (ix < 0 || ix >= xres) continue;
          double rx2 = pow(fabs(fx - (ix + 0.5)) * inv_hpix, 2);
          double u2 = rx2 + ry2;
          if (u2 >= u_max_sq) continue;
          double w = projected_kernel_eval(kt, sqrt(u2));
          int pix = row + ix;
          gas_img[pix] += p->mass * w;
          gtmp_img[pix] += p->mass * p->cooling_data.subgrid_temp * w;
        }
      }
    }

    /* === Stars === */
    for (int j = 0; j < c->stars.count; j++) {
      struct spart *sp = &c->stars.parts[j];

      double P[3] = {sp->x[0] - s->dim[0] * 0.5, sp->x[1] - s->dim[1] * 0.5,
                     sp->x[2] - s->dim[2] * 0.5};
      double V[3] = {P[0] - cam_rel[0], P[1] - cam_rel[1], P[2] - cam_rel[2]};
      double x_cam = V[0] * Rvec[0] + V[1] * Rvec[1] + V[2] * Rvec[2];
      double y_cam = V[0] * Uvec[0] + V[1] * Uvec[1] + V[2] * Uvec[2];
      double z_cam = V[0] * F[0] + V[1] * F[1] + V[2] * F[2];
      if (z_cam <= 1e-6) continue;
      double fx = (x_cam / z_cam / max_x * 0.5 + 0.5) * xres;
      double fy = (y_cam / z_cam / max_y * 0.5 + 0.5) * yres;

      double h_pix = (sp->h / R) * px_per_rad_x;
      double inv_hpix = 1.0 / h_pix;
      if (h_pix < 0.5) {
        int ix = (int)floor(fx + 0.5), iy = (int)floor(fy + 0.5);
        if (ix >= 0 && ix < xres && iy >= 0 && iy < yres)
          star_img[iy * xres + ix] += sp->mass;
        continue;
      }

      int fx_i = (int)floor(fx);
      int fy_i = (int)floor(fy);
      int delta = (int)ceil(u_max * h_pix);

      for (int dj = -delta; dj <= delta; dj++) {
        int iy = fy_i + dj;
        if (iy < 0 || iy >= yres) continue;
        int row = iy * xres;
        double ry2 = pow(fabs(fy - (iy + 0.5)) * inv_hpix, 2);

        for (int di = -delta; di <= delta; di++) {
          int ix = fx_i + di;
          if (ix < 0 || ix >= xres) continue;
          double rx2 = pow(fabs(fx - (ix + 0.5)) * inv_hpix, 2);
          double u2 = rx2 + ry2;
          if (u2 >= u_max_sq) continue;
          double w = projected_kernel_eval(kt, sqrt(u2));
          star_img[row + ix] += sp->mass * w;
        }
      }
    }
  } /* end cells */
}

/**
 * @brief Run the imaging loop for all rotation frames, writing each
 *        rotated view into one HDF5 file using 3-D datasets
 *        of shape (nr_frames, xres, yres), chunked per-frame and
 *        compressed with deflate.
 */
void imaging_compute_angular_images(struct space *s) {
  ticks tic = getticks();

  struct engine *e = s->e;
  struct image_common_data *id = e->image_data;
  const int nf = id->nr_rotation_frames;
  const int xres = id->xres;
  const int yres = id->yres;
  const size_t npix = (size_t)xres * yres;

  /* 3D dataset dims: [frames, X, Y] */
  hsize_t dims3[3] = {(hsize_t)nf, (hsize_t)xres, (hsize_t)yres};

  /* Build temporary & final names */
  char filename[256];
  char final_filename[256];
  snprintf(filename, sizeof(filename), "%s_tmp/%s_%d.hdf5", id->output_dir,
           id->base_name, id->frame_number);
  snprintf(final_filename, sizeof(final_filename), "%s/%s_%d.hdf5",
           id->output_dir, id->base_name, id->frame_number);

  /* 1) Create the HDF5 file */
  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) {
    error("Failed to create HDF5 file %s.", filename);
    return;
  }

  /* 2) Create a property list for chunking+compression */
  hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
  /* One chunk per frame */
  hsize_t chunk3[3] = {1, (hsize_t)xres, (hsize_t)yres};
  H5Pset_chunk(dcpl, 3, chunk3);
  H5Pset_deflate(dcpl, 4);

  /* 3) Create one 3D dataspace */
  hid_t space3 = H5Screate_simple(3, dims3, NULL);

  /* 4) Create four 3D datasets at the file root */
  hid_t d_dm = H5Dcreate2(file_id, "dark_matter", H5T_NATIVE_DOUBLE, space3,
                          H5P_DEFAULT, dcpl, H5P_DEFAULT);
  hid_t d_gas = H5Dcreate2(file_id, "gas", H5T_NATIVE_DOUBLE, space3,
                           H5P_DEFAULT, dcpl, H5P_DEFAULT);
  hid_t d_str = H5Dcreate2(file_id, "stars", H5T_NATIVE_DOUBLE, space3,
                           H5P_DEFAULT, dcpl, H5P_DEFAULT);
  hid_t d_gtmp = H5Dcreate2(file_id, "gas_temperature", H5T_NATIVE_DOUBLE,
                            space3, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  /* We can now close the dataspace & DCPL handles */
  H5Sclose(space3);
  H5Pclose(dcpl);

  /* Prepare a memory‐space for writing one 2D slice [1,xres,yres] */
  hsize_t memdims[3] = {1, (hsize_t)xres, (hsize_t)yres};
  hid_t memspace = H5Screate_simple(3, memdims, NULL);

  /* Remember the base phi just once */
  double phi0 = id->sphere_camera_position[2];

  /* 5) Loop over each rotation frame and write into the 3D datasets */
  for (int f = 0; f < nf; f++) {
    /* Update the spherical camera position for this frame */
    id->sphere_camera_position[2] = phi0 + (2.0 * M_PI * f) / (double)nf;

    /* Allocate & zero each thread’s buffers */
    imaging_allocate_threadimages(e);

    /* Project all cells */
    threadpool_map(&e->threadpool, imaging_cell_mapper, NULL, s->nr_cells, 1,
                   threadpool_auto_chunk_size, e);

    /* Reduce into thread-0 arrays */
    double *dm = id->dm_images[0];
    double *gas = id->gas_images[0];
    double *str = id->star_images[0];
    double *gtmp = id->gas_temp_images[0];
    for (int t = 1; t < e->nr_threads; t++) {
      double *dm_t = id->dm_images[t];
      double *gas_t = id->gas_images[t];
      double *str_t = id->star_images[t];
      double *gtmp_t = id->gas_temp_images[t];
      for (size_t i = 0; i < npix; i++) {
        dm[i] += dm_t[i];
        gas[i] += gas_t[i];
        str[i] += str_t[i];
        gtmp[i] += gtmp_t[i];
      }
    }
    /* Normalize gas temperature */
    for (size_t i = 0; i < npix; i++) {
      if (gas[i] > 0.0)
        gtmp[i] /= gas[i];
      else
        gtmp[i] = 0.0;
    }

    /* Select hyperslab offset = [f,0,0], count = [1,xres,yres] */
    hsize_t offset[3] = {(hsize_t)f, 0, 0};
    hsize_t count[3] = {1, (hsize_t)xres, (hsize_t)yres};

    /* Write each slice */
    {
      hid_t fs = H5Dget_space(d_dm);
      H5Sselect_hyperslab(fs, H5S_SELECT_SET, offset, NULL, count, NULL);
      H5Dwrite(d_dm, H5T_NATIVE_DOUBLE, memspace, fs, H5P_DEFAULT, dm);
      H5Sclose(fs);
    }
    {
      hid_t fs = H5Dget_space(d_gas);
      H5Sselect_hyperslab(fs, H5S_SELECT_SET, offset, NULL, count, NULL);
      H5Dwrite(d_gas, H5T_NATIVE_DOUBLE, memspace, fs, H5P_DEFAULT, gas);
      H5Sclose(fs);
    }
    {
      hid_t fs = H5Dget_space(d_str);
      H5Sselect_hyperslab(fs, H5S_SELECT_SET, offset, NULL, count, NULL);
      H5Dwrite(d_str, H5T_NATIVE_DOUBLE, memspace, fs, H5P_DEFAULT, str);
      H5Sclose(fs);
    }
    {
      hid_t fs = H5Dget_space(d_gtmp);
      H5Sselect_hyperslab(fs, H5S_SELECT_SET, offset, NULL, count, NULL);
      H5Dwrite(d_gtmp, H5T_NATIVE_DOUBLE, memspace, fs, H5P_DEFAULT, gtmp);
      H5Sclose(fs);
    }
  }

  /* Clean up HDF5 handles */
  H5Sclose(memspace);
  H5Dclose(d_dm);
  H5Dclose(d_gas);
  H5Dclose(d_str);
  H5Dclose(d_gtmp);
  H5Fclose(file_id);

  /* Move the temporary file to the final name */
  if (rename(filename, final_filename) != 0) {
    error("Failed to rename %s to %s.", filename, final_filename);
  } else {
    if (e->verbose) {
      message("Wrote %d rotation frames to %s.", nf, final_filename);
    }
  }

  /* Advance frame counter */
  id->frame_number++;

  if (e->verbose) {
    message("Computed %d rotation frames in %.3f %s.", nf,
            clocks_from_ticks(getticks() - tic), clocks_getunit());
  }
}
