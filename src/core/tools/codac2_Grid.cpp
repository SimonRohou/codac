/** 
 *  codac2_Grid.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_serialization.h"
#include "codac2_Grid.h"

using namespace std;
using namespace codac2;

Grid::Grid()
{ }

Grid::Grid(std::ifstream& f)
{
  deserialize(f);
}

Grid::Grid(const Vector& dx, const Vector& origin)
  : Grid(dx, origin, Eigen::Matrix<int,-1,1>::Zero(dx.size()))
{ }

Grid::Grid(const Vector& dx, const Vector& origin, const Eigen::Matrix<int,-1,1>& local_shift)
  : _dx(dx), _origin(origin), _local_shift(local_shift)
{
  assert_release(all_same_size(_dx, _origin, _local_shift));
}

size_t Grid::dim() const
{
  return _dx.size();
}

const Vector& Grid::dx() const
{
  return _dx;
}

const Vector& Grid::origin() const
{
  return _origin;
}

const Eigen::Matrix<int,-1,1>& Grid::local_shift() const
{
  return _local_shift;
}

std::pair<GridBox,GridBox> Grid::bisect(const GridBox& x) const
{
  assert_release(x.size() == (Index)dim());
  assert(x.has_integer_bounds());
  
  auto y = x.bisect_largest();

  assert([&]() {
    Index i = x.max_diam_index();
    return std::floor(y.first[i].ub())+1 == std::ceil(y.second[i].lb());
  } ());

  GridBox b1(y.first.lb(),floor(y.first.ub()).eval());
  GridBox b2(ceil(y.second.lb()).eval(),y.second.ub());

  assert(b1.has_integer_bounds() && b2.has_integer_bounds());
  return { b1, b2 };
}

IntervalVector Grid::to_reals(const GridBox& x, bool local_indices) const
{
  assert(x.size() == (Index)dim());

  if(x.is_empty())
    return IntervalVector::empty(dim());

  assert(x.has_integer_bounds());

  Vector i_lb = x.lb();
  if(local_indices)
    i_lb += _local_shift.template cast<double>();

  // Note: it is important to compute the bounds independently,
  // without any cast to intervals, in order to avoid outward rounding 
  // and so inconsistent grid cells.
  auto w_lb = _origin + _dx.cwiseProduct(i_lb);

  IntervalVector w(w_lb);

  Vector i_ub = x.ub() + Vector::ones(dim());
  if(local_indices)
    i_ub += _local_shift.template cast<double>();

  auto w_ub = _origin + _dx.cwiseProduct(i_ub);

  assert(w_lb != w_ub);
  return { w_lb.eval(), w_ub.eval() };
}

GridBox Grid::to_grid(const IntervalVector& x) const
{
  const Index n = dim();
  assert(x.size() == n); 

  IntervalVector x_ = x - _origin.template cast<Interval>();
  GridBox g(n);

  for(Index i = 0 ; i < n ; i++)
  {
    double lb = (codac2::floor(Interval(x_[i].lb())/_dx[i])).lb();
    g[i] = Interval(lb, std::max(lb,(codac2::ceil(Interval(x_[i].ub())/_dx[i])-1).ub())) - _local_shift[i];
  }

  assert(!g.is_empty());
  assert(g.has_integer_bounds());
  return g;
}

void Grid::serialize(std::ofstream& f) const
{
  assert(f.is_open());
  codac2::serialize(f, _dx);
  codac2::serialize(f, _origin);
  codac2::serialize(f, _local_shift);
}

void Grid::deserialize(std::ifstream& f)
{
  assert(f.is_open());
  codac2::deserialize(f, _dx);
  codac2::deserialize(f, _origin);
  codac2::deserialize(f, _local_shift);
  assert(_dx.min_coeff() > 0);
}