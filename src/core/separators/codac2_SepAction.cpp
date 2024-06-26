/** 
 *  codac2_SepAction.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou, Luc Jaulin
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_SepAction.h"

using namespace std;
using namespace codac2;

std::shared_ptr<Sep> SepAction::copy() const
{
  return std::make_shared<SepAction>(*this);
}

BoxPair SepAction::separate(const IntervalVector& x) const
{
  assert(x.size() == this->size());
  IntervalVector _x(_s(x));
  auto x_sep = _sep.front().separate(_x);
  return { __s(x_sep.inner), __s(x_sep.outer) };
}