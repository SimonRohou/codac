/** 
 *  \file codac2_analytic_variables.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <iostream>
#include "codac2_AnalyticExpr.h"
#include "codac2_VarBase.h"

namespace codac2
{
  template<typename T>
  class AnalyticVarExpr : public AnalyticExpr<T>, public VarBase
  {
    public:

      AnalyticVarExpr()
      { }

      virtual const ExprID& unique_id() const
      {
        return AnalyticExpr<T>::unique_id();
      }

      T fwd_eval(ValuesMap& v, [[maybe_unused]] Index total_input_size) const
      {
        return AnalyticExpr<T>::value(v);
      }
      
      void bwd_eval([[maybe_unused]] ValuesMap& v) const
      { }

      void replace_expr([[maybe_unused]] const ExprID& old_expr_id, [[maybe_unused]] const std::shared_ptr<ExprBase>& new_expr)
      { }
      
      operator std::shared_ptr<AnalyticExpr<T>>() const
      {
        return std::dynamic_pointer_cast<AnalyticExpr<T>>(this->copy());
      }

      virtual bool belongs_to_args_list(const FunctionArgsList& args) const
      {
        for(const auto& xi : args)
          if(xi->unique_id() == this->unique_id())
            return true;
        return false;
      }
  };

  class ScalarVar : public AnalyticVarExpr<ScalarOpValue>
  {
    public:

      ScalarVar()
      { }

      ScalarVar(const ScalarVar& x)
        : AnalyticVarExpr<ScalarOpValue>(x)
      { }

      std::shared_ptr<VarBase> arg_copy() const
      {
        return std::make_shared<ScalarVar>(*this);
      }

      std::shared_ptr<ExprBase> copy() const
      {
        return std::make_shared<ScalarVar>(*this);
      }

      Index size() const
      {
        return 1;
      }

      std::shared_ptr<AnalyticExpr<ScalarOpValue>> operator-() const
      {
        return std::make_shared<AnalyticOperationExpr<SubOp,ScalarOpValue,ScalarOpValue>>(*this);
      }
  };

  class VectorVar : public AnalyticVarExpr<VectorOpValue>
  {
    public:

      explicit VectorVar(Index n)
        : _n(n)
      {
        assert_release(n > 0);
      }

      VectorVar(const VectorVar& x)
        : AnalyticVarExpr<VectorOpValue>(x), _n(x._n)
      { }

      std::shared_ptr<VarBase> arg_copy() const
      {
        return std::make_shared<VectorVar>(*this);
      }

      std::shared_ptr<ExprBase> copy() const
      {
        return std::make_shared<VectorVar>(*this);
      }

      Index size() const
      {
        return _n;
      }

      std::shared_ptr<AnalyticExpr<ScalarOpValue>> operator[](Index i) const
      {
        assert_release(i >= 0 && i < _n);
        return std::make_shared<AnalyticOperationExpr<ComponentOp,ScalarOpValue,VectorOpValue>>(
          std::dynamic_pointer_cast<AnalyticExpr<VectorOpValue>>(this->copy()), i);
      }

      std::shared_ptr<AnalyticExpr<VectorOpValue>> subvector(Index i, Index j) const
      {
        assert_release(i >= 0 && i < _n && j >= i && j < _n);
        return std::make_shared<AnalyticOperationExpr<SubvectorOp,VectorOpValue,VectorOpValue>>(
          std::dynamic_pointer_cast<AnalyticExpr<VectorOpValue>>(this->copy()), i, j);
      }

    protected:

      Index _n;
  };
}