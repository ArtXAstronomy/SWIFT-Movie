/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_REMIX_HYDRO_IO_H
#define SWIFT_REMIX_HYDRO_IO_H

/**
 * @file REMIX/hydro_io.h
 * @brief REMIX implementation of SPH (Sandnes et al. 2025) i/o routines
 */

#include "adiabatic_index.h"
#include "hydro.h"
#include "hydro_parameters.h"
#include "io_properties.h"
#include "kernel_hydro.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void hydro_read_particles(struct part* parts,
                                        struct io_props* list,
                                        int* num_fields) {

#ifdef PLANETARY_FIXED_ENTROPY
  *num_fields = 10;
#else
  *num_fields = 9;
#endif

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                parts, mass);
  list[3] = io_make_input_field("SmoothingLength", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, h);
  list[4] = io_make_input_field("InternalEnergy", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts, u);
  list[5] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_input_field("Accelerations", FLOAT, 3, OPTIONAL,
                                UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] = io_make_input_field("Density", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_DENSITY, parts, rho);
  list[8] = io_make_input_field("MaterialIDs", INT, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, mat_id);
#ifdef PLANETARY_FIXED_ENTROPY
  list[9] = io_make_input_field("Entropies", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_PHYSICAL_ENTROPY_PER_UNIT_MASS, parts,
                                s_fixed);
#endif
}

INLINE static void convert_S(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {

  ret[0] = hydro_get_comoving_entropy(p, xp);
}

INLINE static void convert_P(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {

  ret[0] = hydro_get_comoving_pressure(p);
}

INLINE static void convert_part_pos(const struct engine* e,
                                    const struct part* p,
                                    const struct xpart* xp, double* ret) {
  const struct space* s = e->s;
  if (s->periodic) {
    ret[0] = box_wrap(p->x[0], 0.0, s->dim[0]);
    ret[1] = box_wrap(p->x[1], 0.0, s->dim[1]);
    ret[2] = box_wrap(p->x[2], 0.0, s->dim[2]);
  } else {
    ret[0] = p->x[0];
    ret[1] = p->x[1];
    ret[2] = p->x[2];
  }
  if (e->snapshot_use_delta_from_edge) {
    ret[0] = min(ret[0], s->dim[0] - e->snapshot_delta_from_edge);
    ret[1] = min(ret[1], s->dim[1] - e->snapshot_delta_from_edge);
    ret[2] = min(ret[2], s->dim[2] - e->snapshot_delta_from_edge);
  }
}

INLINE static void convert_part_vel(const struct engine* e,
                                    const struct part* p,
                                    const struct xpart* xp, float* ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology* cosmo = e->cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;
  const float dt_kick_grav_mesh = e->dt_kick_grav_mesh_for_io;

  const integertime_t ti_beg = get_integer_time_begin(ti_current, p->time_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, p->time_bin);

  /* Get time-step since the last kick */
  float dt_kick_grav, dt_kick_hydro;
  if (with_cosmology) {
    dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_current);
    dt_kick_grav -=
        cosmology_get_grav_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
    dt_kick_hydro = cosmology_get_hydro_kick_factor(cosmo, ti_beg, ti_current);
    dt_kick_hydro -=
        cosmology_get_hydro_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
  } else {
    dt_kick_grav = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
    dt_kick_hydro = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
  }

  /* Extrapolate the velocites to the current time (hydro term)*/
  ret[0] = xp->v_full[0] + p->a_hydro[0] * dt_kick_hydro;
  ret[1] = xp->v_full[1] + p->a_hydro[1] * dt_kick_hydro;
  ret[2] = xp->v_full[2] + p->a_hydro[2] * dt_kick_hydro;

  /* Add the gravity term */
  if (p->gpart != NULL) {
    ret[0] += p->gpart->a_grav[0] * dt_kick_grav;
    ret[1] += p->gpart->a_grav[1] * dt_kick_grav;
    ret[2] += p->gpart->a_grav[2] * dt_kick_grav;
  }

  /* And the mesh gravity term */
  if (p->gpart != NULL) {
    ret[0] += p->gpart->a_grav_mesh[0] * dt_kick_grav_mesh;
    ret[1] += p->gpart->a_grav_mesh[1] * dt_kick_grav_mesh;
    ret[2] += p->gpart->a_grav_mesh[2] * dt_kick_grav_mesh;
  }

  /* Conversion from internal units to peculiar velocities */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}

INLINE static void convert_part_potential(const struct engine* e,
                                          const struct part* p,
                                          const struct xpart* xp, float* ret) {

  if (p->gpart != NULL)
    ret[0] = gravity_get_comoving_potential(p->gpart);
  else
    ret[0] = 0.f;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void hydro_write_particles(const struct part* parts,
                                         const struct xpart* xparts,
                                         struct io_props* list,
                                         int* num_fields) {

  *num_fields = 11;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_part(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, parts, xparts,
      convert_part_pos, "Positions of the particles");
  list[1] = io_make_output_field_convert_part(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, parts, xparts,
      convert_part_vel, "Velocities of the particles");
  list[2] = io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts,
                                 mass, "Masses of the particles");
  list[3] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, parts, h,
      "Smoothing lengths (FWHM of the kernel) of the particles");
  list[4] = io_make_output_field(
      "InternalEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      -3.f * hydro_gamma_minus_one, parts, u,
      "Thermal energies per unit mass of the particles");
  list[5] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           parts, id, "Unique IDs of the particles");
  list[6] = io_make_output_field("Densities", FLOAT, 1, UNIT_CONV_DENSITY, -3.f,
                                 parts, rho, "Densities of the particles");
  list[7] = io_make_output_field_convert_part(
      "Entropies", FLOAT, 1, UNIT_CONV_ENTROPY_PER_UNIT_MASS, 0.f, parts,
      xparts, convert_S, "Entropies per unit mass of the particles");
  list[8] =
      io_make_output_field("MaterialIDs", INT, 1, UNIT_CONV_NO_UNITS, 0.f,
                           parts, mat_id, "Material IDs of the particles");
  list[9] = io_make_output_field_convert_part(
      "Pressures", FLOAT, 1, UNIT_CONV_PRESSURE, -3.f * hydro_gamma, parts,
      xparts, convert_P, "Pressures of the particles");
  list[10] = io_make_output_field_convert_part(
      "Potentials", FLOAT, 1, UNIT_CONV_POTENTIAL, 0.f, parts, xparts,
      convert_part_potential, "Gravitational potentials of the particles");
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
INLINE static void hydro_write_flavour(hid_t h_grpsph) {

  io_write_attribute_s(h_grpsph, "Density estimate", "Evolved in time");
  io_write_attribute_s(h_grpsph, "EoM free functions", "Evolved densities");
  io_write_attribute_s(h_grpsph, "Kernel gradients",
                       "Linear-order reproducing kernels with grad-h terms");
  io_write_attribute_s(h_grpsph, "Vacuum boundary treatment",
                       "As presented in Sandnes et al. (2025)");
  io_write_attribute_s(
      h_grpsph, "Artificial viscosity",
      "With linear reconstruction, van Leer slope limiter, Balsara switch");
  io_write_attribute_s(
      h_grpsph, "Artificial diffusion",
      "With linear reconstruction, van Leer slope limiter, Balsara switch");
  io_write_attribute_s(h_grpsph, "Normalising term",
                       "As presented in Sandnes et al. (2025)");
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
INLINE static int writeEntropyFlag(void) { return 0; }

#endif /* SWIFT_REMIX_HYDRO_IO_H */
