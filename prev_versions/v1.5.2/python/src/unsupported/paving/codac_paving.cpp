//============================================================================
//                               P Y I B E X
// File        : pyIbex_GeoImages.cpp
// Author      : Benoit Desrochers
// Copyright   : Benoit Desrochers
// License     : See the LICENSE file
// Created     : Aug 18, 2016
//============================================================================



#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>

#include <codac_IntervalVector.h>

#include <codac_SepPaving.h>
#include <codac_PavingVisitor.h>
#include "../../core/contractors/static/codac_py_Ctc.h"
#include "../../core/separators/codac_py_Sep.h"

#include <codac_Ctc.h>
#include <ibex_Sep.h>



namespace py = pybind11;
using py::self;
using py::init;
using py::keep_alive;
using py::class_;


using namespace pybind11::literals;

using namespace codac;


// using SepVisitor = PavingVisitor<SepPaving, IntervalVector>;
class pySepVisitor : public SepVisitor {
public:
  /* Inherit the constructors */
  using SepVisitor::SepVisitor;
  /* Trampoline (need one for each virtual function) */
  virtual void visit_leaf(const IntervalVector& boxIn, IntervalVector boxOut){
    PYBIND11_OVERLOAD_PURE(void, SepVisitor, visit_leaf, boxIn, boxOut);
  }

  virtual void visit_node(const IntervalVector& boxIn, const IntervalVector& boxOut){
    PYBIND11_OVERLOAD(void, SepVisitor, visit_node, boxIn, boxOut);
  }

  virtual void pre_visit(const SepPaving& p){
    PYBIND11_OVERLOAD(void, SepVisitor, pre_visit, p);
  }
  virtual void post_visit(const SepPaving& p){
    PYBIND11_OVERLOAD(void, SepVisitor, post_visit, p);
  }
};


void export_unsupported_paving(py::module& m, py::class_<Ctc, pyCtc>& ctc, py::class_<Sep, pySep>& sep)
{
  // py::module m("image", "python binding CtcIMage Class");

  class_<SepPaving>(m, "SepPaving", sep)
      // .def(init<int>())
      .def(init<IntervalVector&> (), "X0"_a)
			.def(init<IntervalVector&, LargestFirst& > (), "X0"_a, "bsc"_a)
			.def(init<IntervalVector&, Sep&,double, bool> () , "X0"_a, "sep"_a, "eps"_a, "display"_a=false )
      .def(init<const std::string&> (), "filename"_a)
      .def( "Reunite", &SepPaving::Reunite )
      .def( "Sivia", ( SepPaving& (SepPaving::*) (Sep& ,double)) &SepPaving::Sivia, py::return_value_policy::take_ownership)
      .def( "visit", &SepPaving::visit, py::arg("visitor"))
      .def( "save", &SepPaving::save)
      .def( "getBoundingBox", &SepPaving::getBoundingBox)
      ;

  class_<SepVisitor, std::unique_ptr<SepVisitor>, pySepVisitor>(m, "SepVisitor")
      .def(init<>())
      .def("visit_leaf", &SepVisitor::visit_leaf)
      .def("visit_node", ( void (SepVisitor::*) ( const IntervalVector&, const IntervalVector&) )&SepVisitor::visit_node)
      .def("pre_visit",  &SepVisitor::pre_visit)
      .def("post_visit", &SepVisitor::post_visit)
      ;

  class_<CtcPaving>(m, "CtcPaving", ctc)
      .def(init<SepPaving&>())
      .def("contract", &CtcPaving::contract)
      ;

    // return m.ptr();
}
