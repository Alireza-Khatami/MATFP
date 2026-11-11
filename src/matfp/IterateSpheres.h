// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include <cmath>
#include <utility>
#include <vector>

#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {
// Optional record of one iterate_sphere call (nullptr = off, MATFP as is).
// One record per stage; planes / concave lines of record i are
// [plane_ptr[i], plane_ptr[i+1]) / [cc_ptr[i], cc_ptr[i+1]).
struct IterTrace {
  // INPUT: as given; SNAP: after update_tan_planes_by_sf_mesh;
  // PLANES: after prepare_next_iterate adds q / removes planes (num_itr > 0);
  // SOLVE: after update_sphere_given_plN_opt (planes = the solve's pN, nN);
  // TANGENT: after update_tangent_points_on_tan_pls / _on_cc_lines
  enum Stage { ST_INPUT = 0, ST_SNAP = 1, ST_PLANES = 2, ST_SOLVE = 3, ST_TANGENT = 4 };
  enum PlaneStatus {
    KEPT = 0,
    ADDED_Q = 1,           // nearest point q added in prepare_next_iterate
    REMOVED_SAME_N_Q = 2,  // removed there: same normal as q
    DELETED_SAME_FID = 3,  // removed in tangent update: face already taken
    DELETED_SIMILAR = 4    // removed in tangent update: is_same_tan_pl
  };
  enum Stop {
    RUNNING = -1,
    GOOD = 0,                  // energy / r^2 < break_threshold
    ENERGY_STALL = 1,          // |change of energy / r^2| <= SCALAR_ZERO_N6
    ITER_LIMIT = 2,            // num_itr > iteration_limit
    SMALL_SQ_RADIUS = 3,       // sq_radius < SCALAR_ZERO_N6 (in the check)
    SOLVE_SIZE_MISMATCH = 4,   // pN / nN sizes differ
    SOLVE_NOT_INVERTIBLE = 5,  // 4x4 not invertible
    SOLVE_SMALL_RADIUS = 6,    // solved radius < SCALAR_ZERO_N3
    SOLVE_SMALL_SQ_RADIUS = 7  // solved sq_radius < SCALAR_ZERO_N6
  };
  std::vector<int> stage, iter, q_skipped;
  std::vector<Vector3> pos;
  std::vector<double> sq_radius, energy_over_sq_radius;  // NaN if not computed
  std::vector<int> plane_ptr{0}, cc_ptr{0};
  std::vector<Vector3> plane_point, plane_normal;
  std::vector<int> plane_fid, plane_status;
  std::vector<double> plane_energy;
  std::vector<Vector3> cc_point, cc_normal;
  std::vector<double> cc_energy;
  int stop = RUNNING, num_itr = 0;
  bool is_good = false;

  // mat_p.tan_planes get status KEPT (first_status for plane 0), then removed
  void add(int st, int it, const MVertex& mat_p, double e_over_sq_r = NAN,
           int q_skip = 0, int first_status = KEPT,
           const std::vector<std::pair<TangentPlane, int>>* removed = nullptr);
};

// Default parameters:
// alpha1 = 0.01;  // energy of distance to tangent point
// alpha2 = 1;     // energy of distance to tangent plane
// alpha3 = 1;     // energy of distance to concave line
bool iterate_sphere(const GEO::Mesh& sf_mesh,
                    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
                    const AABBWrapper& aabb_wrapper, MVertex& mat_p,
                    bool is_debug = false, double alpha1 = 0.01,
                    double alpha2 = 1, double alpha3 = 1,
                    const double break_threshold = SCALAR_ZERO_N4,
                    const int iteration_limit = 30,
                    IterTrace* trace = nullptr);
}  // namespace matfp
