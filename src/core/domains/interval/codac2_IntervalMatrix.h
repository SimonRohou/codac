/** 
 *  \file codac2_IntervalMatrix.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2023 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_matrices.h"
#include "codac2_Vector.h"

namespace codac2
{
  using IntervalMatrix = Eigen::Matrix<Interval,-1,-1>;

  inline std::ostream& operator<<(std::ostream& os, const IntervalMatrix& x)
  {
    if(x.is_empty())
      return os << "( empty matrix )";

    else
      return operator<<(os,
        static_cast<const Eigen::DenseBase<IntervalMatrix>&>(x));

    #if 0 // IBEX style
      os << "(";
      for(Index i = 0 ; i < x.rows() ; i++)
      {
        os << (i!=0 ? " " : "") << "(";
        for(Index j = 0 ; j < x.cols() ; j++)
          os << x(i,j) << (j<x.cols()-1 ? " ; " : "");
        os << ")";
        if(i < x.rows()-1) os << std::endl;
      }
      os << ")";
      return os;
    #endif
  }
}