/** 
 *  \file codac2_CtcUnion.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <type_traits>
#include "codac2_CtcWrapper.h"
#include "codac2_Collection.h"
#include "codac2_template_tools.h"

namespace codac2
{
  template<typename X>
  class CtcUnion : public Ctc<CtcUnion<X>,X>
  {
    public:

      explicit CtcUnion(size_t n)
        : Ctc<CtcUnion<X>,X>(n)
      {
        if constexpr(std::is_same_v<X,Interval>)
        {
          assert(n == 1);
        }
      }

      template<typename C, typename = typename std::enable_if<(
          (std::is_base_of_v<CtcBase<X>,C> && !std::is_same_v<CtcUnion<X>,C>) || std::is_same_v<std::shared_ptr<CtcBase<X>>,C>
        ), void>::type>
      CtcUnion(const C& c)
        : Ctc<CtcUnion<X>,X>(size_of(c)), _ctcs(c)
      { }

      template<typename... C, typename = typename std::enable_if<(true && ... && (
          (std::is_base_of_v<CtcBase<X>,C> || std::is_same_v<std::shared_ptr<CtcBase<X>>,C>)
        )), void>::type>
      CtcUnion(const C&... c)
        : Ctc<CtcUnion<X>,X>(size_first_item(c...)), _ctcs(c...)
      {
        assert_release(all_same_size(c...));
      }

      void contract(X& x) const
      {
        auto result = x;
        result.set_empty();

        for(const auto& ci : _ctcs)
        {
          auto saved_x = x;
          ci->contract(saved_x);
          result |= saved_x;
        }

        x = result;
      }

      template<typename C, typename = typename std::enable_if<
          std::is_base_of_v<CtcBase<X>,C>
        >::type>
      CtcUnion<X>& operator|=(const C& c)
      {
        assert_release(c.size() == this->size());
        _ctcs.add_shared_ptr(std::make_shared<C>(c));
        return *this;
      }

      CtcUnion<X>& operator|=(const std::shared_ptr<CtcBase<X>>& c)
      {
        assert_release(c->size() == this->size());
        _ctcs.add_shared_ptr(c);
        return *this;
      }

    protected:

      Collection<CtcBase<X>> _ctcs;
  };

  template<typename C1, typename C2, typename = typename std::enable_if<(
      std::is_base_of_v<CtcBase<typename C1::ContractedType>,C1> &&
      std::is_base_of_v<CtcBase<typename C1::ContractedType>,C2> &&
      std::is_same_v<typename C1::ContractedType,typename C2::ContractedType>
    )>>
  inline CtcUnion<typename C1::ContractedType> operator|(const C1& c1, const C2& c2)
  {
    return CtcUnion<typename C1::ContractedType>(c1,c2);
  }

  template<typename C2, typename = typename std::enable_if<
      std::is_base_of_v<CtcBase<IntervalVector>,C2>
    >::type>
  inline CtcUnion<IntervalVector> operator|(const IntervalVector& c1, const C2& c2)
  {
    assert_release(c1.size() == c2.size());
    return CtcUnion<IntervalVector>(CtcWrapper_<IntervalVector>(c1),c2);
  }

  template<typename C1, typename = typename std::enable_if<
      std::is_base_of_v<CtcBase<IntervalVector>,C1>
    >::type>
  inline CtcUnion<IntervalVector> operator|(const C1& c1, const IntervalVector& c2)
  {
    assert_release(c1.size() == c2.size());
    return CtcUnion<IntervalVector>(c1,CtcWrapper_<IntervalVector>(c2));
  }
}