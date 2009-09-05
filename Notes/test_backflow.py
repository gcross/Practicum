#@+leo-ver=4-thin
#@+node:gcross.20090904201537.1554:@thin test_backflow.py
#@+others
#@+node:gcross.20090904201537.1559:Functions
#@+node:gcross.20090904201537.1558:make_symbols
def make_x(number_of_particles,number_of_dimensions):
    return array([[Symbol("x_{0}_{1}".format(i,k)) for k in xrange(number_of_dimensions)] for i in xrange(number_of_particles)])
#@-node:gcross.20090904201537.1558:make_symbols
#@+node:gcross.20090904201537.1790:generate_random_substitutions
def generate_random_substitutions(*symbols):
    substitutions = {}
    for symbol in symbols:
        substitutions[symbol] = random()
    return substitutions
#@-node:gcross.20090904201537.1790:generate_random_substitutions
#@-node:gcross.20090904201537.1559:Functions
#@+node:gcross.20090904201537.1798:Formula Factories
#@+node:gcross.20090904201537.1561:make_r
def make_r(x):
    def r(i,j):
        return sqrt(pysum((x[i,:]-x[j,:])**2))
    return r
#@-node:gcross.20090904201537.1561:make_r
#@+node:gcross.20090904201537.1562:make_rho
def make_rho(x):
    def rho(i):
        return sqrt(pysum(x[i,0:1]**2))
    return rho
#@-node:gcross.20090904201537.1562:make_rho
#@+node:gcross.20090904201537.1794:make_diff_log_trial
def make_diff_log_trial(a,r,x):
    def diff_log_trial(i,k):
        return pysum(a*(x[i,k]-x[j,k])/(r(i,j)**3 * (1-a/r(i,j))) for j in xrange(x.shape[0]) if j != i)
    return diff_log_trial
#@nonl
#@-node:gcross.20090904201537.1794:make_diff_log_trial
#@+node:gcross.20090904201537.1799:backflow
def backflow(a,r,rho,x):
    number_of_particles = x.shape[0]
    return \
        pysum(
            a*(x[i,0]*x[j,1]-x[i,1]*x[j,0])/rho(i)**2/r(i,j)**3/(1-a/r(i,j))
            for i in xrange(number_of_particles) for j in xrange(number_of_particles) if i != j
        )
#@-node:gcross.20090904201537.1799:backflow
#@-node:gcross.20090904201537.1798:Formula Factories
#@-others

#@<< Import needed modules >>
#@+node:gcross.20090904201537.1564:<< Import needed modules >>
import unittest
from paycheck import *
from numpy import array, dot
from sympy import *
import itertools
from random import random, randint
from __builtin__ import sum as pysum
#@-node:gcross.20090904201537.1564:<< Import needed modules >>
#@nl

#@<< Tests >>
#@+node:gcross.20090904201537.1557:<< Tests >>
#@+others
#@+node:gcross.20090904201537.1795:(verified)
#@+node:gcross.20090904201537.1563:test_diff_r
@with_checker
def test_diff_r(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(1,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    r = make_r(x)
    i = randint(0,number_of_particles-2)
    j = randint(i+1,number_of_particles-1)
    if number_of_dimensions == 1:
        k = 0
    else:
        k = randint(0,number_of_dimensions-1)

    substitutions = generate_random_substitutions(*x.ravel())
    self.assertAlmostEqual(r(i,j).diff(x[i,k]).subs(substitutions),(+(x[i,k]-x[j,k])/r(i,j)).subs(substitutions))
    self.assertAlmostEqual(r(i,j).diff(x[j,k]).subs(substitutions),(-(x[i,k]-x[j,k])/r(i,j)).subs(substitutions))
#@-node:gcross.20090904201537.1563:test_diff_r
#@+node:gcross.20090904201537.1789:test_diff_log_trial
@with_checker
def test_diff_log_trial(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(1,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)

    log_trial = 0
    for i in xrange(number_of_particles):
        for j in xrange(i+1,number_of_particles):
            log_trial += log(1-a/r(i,j))

    i = randint(0,number_of_particles-1)
    k = randint(0,number_of_dimensions-1)

    computed_diff = log_trial.diff(x[i,k])
    formula_diff = make_diff_log_trial(a,r,x)(i,k)

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    self.assertAlmostEqual(computed_diff.subs(substitutions),formula_diff.subs(substitutions))
#@-node:gcross.20090904201537.1789:test_diff_log_trial
#@+node:gcross.20090904201537.1797:test_backflow
@with_checker(number_of_calls=100)
def test_backflow(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    diff_log_trial = make_diff_log_trial(a,r,x)

    gradient_log_trial = array([[diff_log_trial(i,k) for k in xrange(2)] for i in xrange(number_of_particles)])

    gradient_feynman = array([[x[i,1]/rho(i)**2,-x[i,0]/rho(i)**2] for i in xrange(number_of_particles)])

    computed_backflow = dot(gradient_log_trial.ravel(),gradient_feynman.ravel())
    formula_backflow = backflow(a,r,rho,x)

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    self.assertAlmostEqual(computed_backflow.subs(substitutions),formula_backflow.subs(substitutions))
#@-node:gcross.20090904201537.1797:test_backflow
#@-node:gcross.20090904201537.1795:(verified)
#@+node:gcross.20090904201537.1560:TestContainer
class TestContainer(unittest.TestCase):
    pass
    #@    @+others
    #@-others
#@-node:gcross.20090904201537.1560:TestContainer
#@-others

tests = [unittest.defaultTestLoader.loadTestsFromTestCase(test_case) for test_case in
    [
        TestContainer,
    ]]
#@-node:gcross.20090904201537.1557:<< Tests >>
#@nl

unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(tests))
#@-node:gcross.20090904201537.1554:@thin test_backflow.py
#@-leo
