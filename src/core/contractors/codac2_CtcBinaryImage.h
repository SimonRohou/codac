/** 
 *  \file codac2_CtcBinaryImage.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou, inspired from the work of Benoit Desrochers, Jan Sliwka
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <iomanip>
#include "codac2_Ctc.h"
#include "codac2_BinaryImage.h"

namespace codac2
{
  template<size_t N=2>
  class CtcBinaryImage : public Ctc<CtcBinaryImage<N>,IntervalVector>
  {
    public:

      CtcBinaryImage(const BinaryImage<N>& img)
        : Ctc<CtcBinaryImage<N>,IntervalVector>(N), _img(img), _real_box(img.real_box())
      { }

      inline BinaryImage<N>& image()
      {
        return _img;
      }

      const IntervalVector& real_box() const
      {
        return _real_box;
      }

      inline void contract(IntervalVector& x) const
      {
        if(!x.intersects(_real_box))
          return;

        else
        {
          if(!_img.data_has_been_loaded())
            _img.load_data();
          _img.convert_as_integrale_img();

          if(x.is_subset(_real_box))
          {
            GridBox i = _img.to_grid(x);
            _img.contract_indices(i);
            x &= _img.to_reals(i);
          }

          else
          {
            // Only the part inside the real_box is contracted
            auto x_img = x & _real_box;
            contract(x_img);
            for(const auto& yi : x.diff(_real_box))
              x_img |= yi;
            x = x_img;
          }
        }
      }

    //protected:

      mutable BinaryImage<N> _img;
      const IntervalVector _real_box;
  };
}