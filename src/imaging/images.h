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
#ifndef SWIFT_IMAGES_H
#define SWIFT_IMAGES_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "engine.h"
#include "lightcone/projected_kernel.h"
#include "parser.h"

/**
 * @brief The data structure for the image common data.
 *
 * This structure contains all the information needed to create an image
 * from the simulation data.
 *
 */
struct image_common_data {

  /*! The number of images */
  int num_images;

  /*! The images themselves */
  struct image_data *images;

  /*! Image resolution along x and y axes */
  int xres;
  int yres;

  /*! The image's pixel size */
  double pixel_size[2];

  /*! Are we doing a slice? */
  int slice;

  /*! Image "thickness" along z axis, by default the whole box will be used */
  double slice_thickness;

  /*! Output dir */
  char output_dir[256];

  /*! Projected kernel lookup table */
  struct projected_kernel_table *projected_kernel_table;

  /*! Are we writing pngs? */
  int write_pngs;

  /*! Are we writing raw arrays? */
  int write_raw_arrays;

  /*! The image lower left corner. */
  double origin[3];

  /*! Are we imaging a subvolume? */
  int subvolume;

  /*! The centre of the subvolume (only applicable if subvolume is set) */
  double subvolume_centre[3];

  /*! The size of the subvolume (only applicable if subvolume is set) */
  double fov[3];
};

/**
 * @brief The data structure for the image data.
 *
 * This structure contains all the information needed to create an image
 * from the simulation data.
 *
 * There can be multiple images in a single simulation, so this structure
 * is generic. The number of images is set by the number of Image blocks
 * and the ImageCommon block in the parameter file.
 */
struct image_data {

  /*! Image resolution along x and y axes */
  int xres;
  int yres;

  /*! Are we doing a slice? */
  int slice;

  /*! Image "thickness" along z axis, by default the whole box will be used */
  double slice_thickness;

  /*! The image's pixel size */
  double pixel_size[2];

  /*! What particle type is this image for? */
  int particle_type;

  /*! What particle field are we using? */
  char field_name[256];

  /*! Output dir */
  char output_dir[256];

  /*! Base name for the image */
  char base_name[256];

  /*! Current frame number */
  int frame_number;

  /*! The index of the image this one is weight by (-1 if unweighted) */
  int weight_by;
};

void imaging_init(struct image_common_data *image_data,
                  struct swift_params *parameter_file, const int verbose,
                  const double dim[3], const int nodeID);
void imaging_write_images(struct engine *e);
int imaging_cell_overlaps_fov(const struct image_common_data *image_data,
                              const struct cell *c);
void imaging_clean(struct image_common_data *image_data);

#endif /* SWIFT_IMAGES_H */
