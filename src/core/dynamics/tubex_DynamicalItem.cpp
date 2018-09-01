/* ============================================================================
 *  tubex-lib - DynamicalItem class
 * ============================================================================
 *  Copyright : Copyright 2017 Simon Rohou
 *  License   : This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 *
 *  Author(s) : Simon Rohou
 *  Bug fixes : -
 *  Created   : 2015
 * ---------------------------------------------------------------------------- */

#include "tubex_DynamicalItem.h"

using namespace std;
using namespace ibex;

namespace tubex
{
  bool DynamicalItem::valid_domain(const ibex::Interval& domain)
  {
    return !domain.is_empty() && !domain.is_unbounded() && !domain.is_degenerated();
  }
}