/* ============================================================================
 *  tubex-lib - Serialization tools
 * ============================================================================
 *  Copyright : Copyright 2017 Simon Rohou
 *  License   : This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 *
 *  Author(s) : Simon Rohou
 *  Bug fixes : -
 *  Created   : 2017
 * ---------------------------------------------------------------------------- */

#ifndef TRAJECTORY_SERIALIZATION_HEADER
#define TRAJECTORY_SERIALIZATION_HEADER

#include <fstream>
#include "tubex.h"

namespace tubex
{
  void serializeTrajectory(std::ofstream& bin_file, const Trajectory& traj, int version_number = SERIALIZATION_VERSION);
  void deserializeTrajectory(std::ifstream& bin_file, Trajectory& traj);
}

#endif