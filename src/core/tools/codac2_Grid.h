/** 
 *  \file codac2_Grid.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <fstream>
#include "codac2_IntervalVector.h"

namespace codac2
{
  using GridBox = IntervalVector;

  class Grid
  {
    public:

      Grid(std::ifstream& f);
      Grid(const Vector& dx, const Vector& origin);
      Grid(const Vector& dx, const Vector& origin, const Eigen::Matrix<int,-1,1>& local_shift);

      size_t dim() const;

      const Vector& dx() const;
      const Vector& origin() const;
      const Eigen::Matrix<int,-1,1>& local_shift() const;

      std::pair<GridBox,GridBox> bisect(const GridBox& x) const;
      
      IntervalVector to_reals(const GridBox& x, bool local_indices = true) const;
      GridBox to_grid(const IntervalVector& x) const;

      void serialize(std::ofstream& f) const;

    protected:

      Grid();
      void deserialize(std::ifstream& f);

    protected:

      Vector _dx; // in real values
      Vector _origin; // in real values
      Eigen::Matrix<int,-1,1> _local_shift; // in indices
  };
}