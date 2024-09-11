from codac import *

# As an example, defining a custom separator as the complementary of another one

class MySep(Sep):

  def __init__(self, S_):
    Sep.__init__(self, S_.size())
    self.S = S_

  def separate(self, x):
    inner,outer = self.S.separate(x) # separates x into two boxes
    return [outer,inner] # inner/outer boxes are permuted

# Trying this custom separator

x = VectorVar(2)
f = AnalyticFunction([x], sqrt(sqr(x[0])+sqr(x[1])))
s = SepInverse(f, [2,4])

a = MySep(s) & MySep(SepWrapper([[-1.5,1.5],[-0.25,0.25]]))
draw_while_paving([[-5,5],[-5,5]], a, 0.1)