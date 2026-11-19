/** 
 *  \file codac2_AnalyticExpr.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou, Damien Massé
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <memory>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>
#include "codac2_ExprBase.h"
#include "codac2_Domain.h"
#include "codac2_FunctionArgsList.h"
#include "codac2_AnalyticType.h"

namespace codac2
{
  // Built once per AnalyticFunction, after its expression and arguments are finalized.
  // ExprID values are global and may be sparse, hence the local compact indexing.
  class AnalyticValuesLayout
  {
    public:
      AnalyticValuesLayout(const FunctionArgsList& args, const std::shared_ptr<ExprBase>& expr)
      {
        for(const auto& arg : args)
          _indices.emplace(arg->unique_id().id(), _indices.size());
        add(expr);
      }

      std::size_t size() const
      {
        return _indices.size();
      }

      std::size_t index_of(const ExprID& id) const
      {
        const auto it = _indices.find(id.id());
        assert(it != _indices.end() && "expression not present in analytic layout");
        return it->second;
      }

    private:
      void add(const std::shared_ptr<ExprBase>& expr)
      {
        if(!expr)
          return;

        const auto [it, inserted] = _indices.emplace(expr->unique_id().id(), _indices.size());
        if(!inserted)
          return;

        for(const auto& child : expr->children_expr_base())
          add(child);
      }

      std::unordered_map<Index,std::size_t> _indices;
  };

  // Each evaluation owns its values. Slots are contiguous and indexed by the
  // immutable layout, rather than individually allocated in a std::map.
  class ValuesMap
  {
    public:
      explicit ValuesMap(const AnalyticValuesLayout& layout)
        : _layout(layout), _values(layout.size())
      { }

      template<typename T,typename U>
      T& set(const ExprID& id, U&& x)
      {
        auto& slot = _values[_layout.index_of(id)];
        if(auto* p = std::get_if<T>(&slot))
        {
          *p = std::forward<U>(x);
          return *p;
        }

        assert(std::holds_alternative<std::monostate>(slot) && "analytic value type mismatch");
        return slot.emplace<T>(std::forward<U>(x));
      }

      template<typename T>
      T& get(const ExprID& id)
      {
        auto* p = std::get_if<T>(&_values[_layout.index_of(id)]);
        assert(p && "missing or mismatched analytic value");
        return *p;
      }

      template<typename T>
      const T& get(const ExprID& id) const
      {
        const auto* p = std::get_if<T>(&_values[_layout.index_of(id)]);
        assert(p && "missing or mismatched analytic value");
        return *p;
      }

    private:
      using Value = std::variant<std::monostate,ScalarType,VectorType,MatrixType>;
      const AnalyticValuesLayout& _layout;
      std::vector<Value> _values;
  };

  template<typename T>
  class AnalyticExpr : public ExprBase
  {
    public:

      virtual const T& fwd_eval(ValuesMap& v, Index total_input_size, bool natural_eval) const = 0;
      virtual void bwd_eval(ValuesMap& v) const = 0;
      virtual std::pair<Index,Index> output_shape() const = 0;

      const T& init_value(ValuesMap& v, const T& x) const
      {
        return v.set<T>(unique_id(), x);
      }

      const T& init_value(ValuesMap& v, T&& x) const
      {
        return v.set<T>(unique_id(), std::move(x));
      }

      T& value(ValuesMap& v) const
      {
        return v.get<T>(unique_id());
      }

      virtual bool belongs_to_args_list(const FunctionArgsList& args) const = 0;
      virtual std::string str(bool in_parentheses = false) const = 0;
      virtual bool is_str_leaf() const = 0;
  };

  template<typename C,typename Y,typename... X>
  class AnalyticOperationExpr : public AnalyticExpr<Y>, public OperationExprBase<AnalyticExpr<X>...>
  {
    public:

      AnalyticOperationExpr(std::shared_ptr<AnalyticExpr<X>>... x)
        : OperationExprBase<AnalyticExpr<X>...>(x...)
      { }

      AnalyticOperationExpr(const AnalyticOperationExpr<C,Y,X...>& e)
        : OperationExprBase<AnalyticExpr<X>...>(e)
      { }

      std::shared_ptr<ExprBase> copy() const override
      {
        return std::make_shared<AnalyticOperationExpr<C,Y,X...>>(*this);
      }

      void replace_arg(const ExprID& old_arg_id, const std::shared_ptr<ExprBase>& new_expr) override
      {
        return OperationExprBase<AnalyticExpr<X>...>::replace_arg(old_arg_id, new_expr);
      }

      const Y& fwd_eval(ValuesMap& v, Index total_input_size, bool natural_eval) const override
      {
        return std::apply(
          [this,&v,total_input_size,natural_eval](auto &&... x) -> const Y&
          {
            if(natural_eval)
              return AnalyticExpr<Y>::init_value(v,
                C::fwd_natural(x->fwd_eval(v, total_input_size, natural_eval)...));

            else
              return AnalyticExpr<Y>::init_value(v,
                C::fwd_centered(x->fwd_eval(v, total_input_size, natural_eval)...));
          },
        this->_x);
      }

      void bwd_eval(ValuesMap& v) const override
      {
        const auto& y = AnalyticExpr<Y>::value(v);

        std::apply([&v,&y](auto &&... x)
        {
          C::bwd(y.a, x->value(v).a...);
        }, this->_x);

        std::apply([&v](auto &&... x)
        {
          (x->bwd_eval(v), ...);
        }, this->_x);
      }

      virtual std::string str(bool in_parentheses = false) const override
      {
        std::string s = std::apply([](auto &&... x) {
          return C::str(x...);
        }, this->_x);
        return in_parentheses ? "(" + s + ")" : s;
      }

      virtual bool is_str_leaf() const override
      {
        return false;
      }

      std::pair<Index,Index> output_shape() const override
      {
        std::pair<Index,Index> s;
        std::apply([&s](auto &&... x)
        {
          s = C::output_shape(x...);
        }, this->_x);
        return s;
      }

      virtual bool belongs_to_args_list(const FunctionArgsList& args) const override
      {
        bool b = true;

        std::apply([&b,&args](auto &&... x)
        {
          ((b &= x->belongs_to_args_list(args)), ...);
        }, this->_x);

        return b;
      }

      std::vector<std::shared_ptr<ExprBase>> children_expr_base() const override
      {
        return OperationExprBase<AnalyticExpr<X>...>::children_expr_base();
      }
  };
}
