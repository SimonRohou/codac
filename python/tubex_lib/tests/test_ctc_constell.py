#!/usr/bin/env python

import unittest
from pyibex import Interval, IntervalVector
from tubex_lib import *
import tubex_lib as tubex

class TestCtcConstell(unittest.TestCase):

  def test_raw(self):

    v_map = [IntervalVector(2)] * 4
    v_map[0] = IntervalVector([1.5,2.5])
    v_map[1] = IntervalVector([2.5,1.5])
    v_map[2] = IntervalVector([4.,2.])
    v_map[3] = IntervalVector([1.,0.5])

    ctc_constell = CtcConstell(v_map)

    x1 = IntervalVector([[0.5,6],[0,3.5]])
    x2 = IntervalVector([[3.5,5.5],[0.5,2.5]])
    x3 = IntervalVector([[2,5],[1,3]])
    x4 = IntervalVector([[1,3],[1,2]])
    x5 = IntervalVector([[1.5,2],[2.5,3]])
    x6 = IntervalVector([[6.5,7.5],[2.5,3.5]])

    ctc_constell.contract(x1)
    ctc_constell.contract(x2)
    ctc_constell.contract(x3)
    ctc_constell.contract(x4)
    ctc_constell.contract(x5)
    ctc_constell.contract(x6)

    self.assertEqual(x1, IntervalVector([[1,4],[0.5,2.5]]))
    self.assertEqual(x2, IntervalVector([4,2]))
    self.assertEqual(x3, IntervalVector([[2.5,4],[1.5,2]]))
    self.assertEqual(x4, IntervalVector([2.5,1.5]))
    self.assertEqual(x5, IntervalVector([1.5,2.5]))
    self.assertEqual(x6, IntervalVector.empty(2))

  def test_in_CN(self):

    v_map = [IntervalVector(2)] * 4
    v_map[0] = IntervalVector([1.5,2.5])
    v_map[1] = IntervalVector([2.5,1.5])
    v_map[2] = IntervalVector([4.,2.])
    v_map[3] = IntervalVector([1.,0.5])

    ctc_constell = CtcConstell(v_map)

    x1 = IntervalVector([[0.5,6],[0,3.5]])
    x2 = IntervalVector([[3.5,5.5],[0.5,2.5]])
    x3 = IntervalVector([[2,5],[1,3]])
    x4 = IntervalVector([[1,3],[1,2]])
    x5 = IntervalVector([[1.5,2],[2.5,3]])
    x6 = IntervalVector([[6.5,7.5],[2.5,3.5]])

    cn = ContractorNetwork()
    cn.add(ctc_constell, [x1])
    cn.add(ctc_constell, [x2])
    cn.add(ctc_constell, [x3])
    cn.add(ctc_constell, [x4])
    cn.add(ctc_constell, [x5])
    cn.add(ctc_constell, [x6])
    cn.contract()

    self.assertEqual(x1, IntervalVector([[1,4],[0.5,2.5]]))
    self.assertEqual(x2, IntervalVector([4,2]))
    self.assertEqual(x3, IntervalVector([[2.5,4],[1.5,2]]))
    self.assertEqual(x4, IntervalVector([2.5,1.5]))
    self.assertEqual(x5, IntervalVector([1.5,2.5]))
    self.assertEqual(x6, IntervalVector.empty(2))

if __name__ ==  '__main__':
  unittest.main()