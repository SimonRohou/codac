/** 
 *  \file
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Simon Rohou, Benoît Desrochers
 *  \copyright  Copyright 2021 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <codac2_Interval.h>
#include <codac2_AnalyticFunction.h>
#include <codac2_analytic_values.h>
#include "codac2_py_AnalyticFunction.h"
#include "codac2_py_CtcInverse.h"
#include "codac2_py_SepInverse.h"
#include "codac2_py_core.h"

using namespace codac2;
namespace py = pybind11;

// 3rd

// contractors
py::class_<Ctc_<IntervalVector>,pyCtcIntervalVector> export_CtcIntervalVector(py::module& m);
void export_CtcInverse(py::module& m, const std::string& export_name, py::class_<Ctc_<IntervalVector>,pyCtcIntervalVector>& ctc);

// domains
py::class_<Interval> export_Interval(py::module& m);
void export_Interval_operations(py::module& m, py::class_<Interval>& py_Interval);
py::class_<IntervalVector> export_IntervalVector(py::module& m);
void export_IntervalVector_operations(py::module& m, py::class_<IntervalVector>& py_IntervalVector);

// functions
void export_ExprWrapperBase(py::module& m);
void export_ScalarVar(py::module& m);
void export_VectorVar(py::module& m);
void export_expression_operations(py::module& m);

// graphics
void export_Figure2D(py::module& m);
void export_StyleProperties(py::module& m);

// paver
void export_Paver(py::module& m, py::class_<Ctc_<IntervalVector>,pyCtcIntervalVector>& ctc);

// separators
py::class_<Sep,pySep> export_Sep(py::module& m);
void export_SepInverse(py::module& m, const std::string& export_name, py::class_<Sep,pySep>& sep);
void export_SepPolygon(py::module& m, py::class_<Sep,pySep>& sep);

// tools
void export_Approx(py::module& m);

// types
py::class_<Vector> export_Vector(py::module& m);


PYBIND11_MODULE(core, m)
{
  m.doc() = string(FOR_MATLAB ? "Matlab" : "Python") + " binding of Codac (core)";
  m.attr("oo") = oo;

  // 3rd

  // contractors
  auto py_ctc = export_CtcIntervalVector(m);
  export_CtcInverse<Interval>(m,"CtcInverse_Interval",py_ctc);
  export_CtcInverse<IntervalVector>(m,"CtcInverse_IntervalVector",py_ctc);

  // domains
  auto py_Interval = export_Interval(m);
  export_Interval_operations(m, py_Interval);
  auto py_IntervalVector = export_IntervalVector(m);
  export_IntervalVector_operations(m, py_IntervalVector);

  // function
  export_ExprWrapperBase(m);
  export_ScalarExpr(m);
  export_VectorExpr(m);
  export_AnalyticFunction<ScalarOpValue>(m,"AnalyticFunction_Scalar");
  export_AnalyticFunction<VectorOpValue>(m,"AnalyticFunction_Vector");
  export_ScalarVar(m);
  export_VectorVar(m);
  export_expression_operations(m);

  // graphics
  export_StyleProperties(m);
  export_Figure2D(m);

  // paver
  export_Paver(m,py_ctc);

  // separators
  auto py_sep = export_Sep(m);
  export_SepInverse<Interval>(m,"SepInverse_Interval",py_sep);
  export_SepInverse<IntervalVector>(m,"SepInverse_IntervalVector",py_sep);
  export_SepPolygon(m,py_sep);

  // tools
  export_Approx(m);

  // types
  auto py_Vector = export_Vector(m);

}