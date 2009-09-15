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
        return sqrt(pysum(x[i,0:2]**2))
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
#@+node:gcross.20090907141528.1358:backflow_numerator
def backflow_numerator(x,i,j):
    number_of_particules = x.shape[0]
    if i == j:
        return 0 * x[0,0]
    else:
        return x[i,0]*x[j,1]-x[i,1]*x[j,0]
#@-node:gcross.20090907141528.1358:backflow_numerator
#@+node:gcross.20090907141528.1360:backflow_denominator
def backflow_denominator(a,r,rho,x,i,j):
    number_of_particules = x.shape[0]
    return rho(i)**2 * r(i,j)**3 * (1-a/r(i,j))
#@-node:gcross.20090907141528.1360:backflow_denominator
#@+node:gcross.20090905201409.1346:backflow_gradient_denominator
def gradient_backflow_denominator(a,r,rho,x,i,j,ip,kp):
    if i == j:
        return 0.0 * x[ip,kp]

    number_of_particles = x.shape[0]
    formula = 0.0 * x[ip,kp]

    r_ij = r(i,j)
    rho_i = rho(i)
    if i==ip and (kp == 0 or kp == 1):
        formula += 2.0 * x[ip,kp] * r_ij**3 * (1.0-a/r_ij)
    if (i==ip):
        factor = +1.0
    elif (j==ip):
        factor = -1.0
    else:
        factor = 0.0
    if factor:
        formula += factor*(rho_i**2 * 3.0 * r_ij * (x[i,kp]-x[j,kp]) * (1.0-a/r_ij) + rho_i**2 * a * (x[i,kp]-x[j,kp])) 

    return formula
#@-node:gcross.20090905201409.1346:backflow_gradient_denominator
#@+node:gcross.20090907141528.1351:backflow_gradient_numerator
def gradient_backflow_numerator(a,r,rho,x,i,j,ip,kp):
    if i == j:
        return 0.0 * x[ip,kp]

    d_kp0 = 1 if kp == 0 else 0
    d_kp1 = 1 if kp == 1 else 0
    d_ipi = 1 if ip == i else 0
    d_ipj = 1 if ip == j else 0

    return (x[j,1]*d_kp0 - x[j,0]*d_kp1)*d_ipi - (x[i,1]*d_kp0 - x[i,0]*d_kp1)*d_ipj
#@-node:gcross.20090907141528.1351:backflow_gradient_numerator
#@+node:gcross.20090908092630.1346:backflow_gradient
def backflow_gradient(a,r,rho,x,ip,kp):

    number_of_particles = x.shape[0]

    def C(i,j):
        return x[i,0]*x[j,1]-x[i,1]*x[j,0]

    def D(i,j):
        return r(i,j)**3 * (1-a/r(i,j))

    formula = 0
    d_kp0 = 1 if kp == 0 else 0
    d_kp1 = 1 if kp == 1 else 0
    for j in xrange(number_of_particles):
        if j == ip:
            continue
        C_ipj = C(ip,j)
        D_ipj = D(ip,j)
        rho_ip = rho(ip)
        rho_j = rho(j)
        r_ipj = r(ip,j)

        formula += 1/D_ipj * (
            (x[j,1]*d_kp0 - x[j,0]*d_kp1)*(1/rho_ip**2 - 1/rho_j**2)
          - C_ipj * ((2*x[ip,kp]*(d_kp0+d_kp1)/rho_ip**4 + (x[ip,kp]-x[j,kp])/r_ipj**2 * (1/rho_ip**2 - 1/rho_j**2) * (3 + a/(r_ipj * (1 - a/r_ipj)))))
        )

    formula *= a

    return formula
#@-node:gcross.20090908092630.1346:backflow_gradient
#@-node:gcross.20090904201537.1798:Formula Factories
#@-others

try:
    import psyco
    psyco.full()
    print "(Loaded psycho)"
except ImportError:
    pass

#@<< Import needed modules >>
#@+node:gcross.20090904201537.1564:<< Import needed modules >>
import unittest
from paycheck import *
from numpy.random import rand
from sympy import *
from numpy import array, dot, zeros, double
import itertools
from random import random, randint
from __builtin__ import sum as pysum

import vpif
#@-node:gcross.20090904201537.1564:<< Import needed modules >>
#@nl

#@<< Tests >>
#@+node:gcross.20090904201537.1557:<< Tests >>
#@+others
#@+node:gcross.20090904201537.1795:(verified)
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
#@+node:gcross.20090906145419.1341:test_diff_r
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
#@-node:gcross.20090906145419.1341:test_diff_r
#@+node:gcross.20090906145419.1345:test_diff_r_cubed
@with_checker
def test_diff_r_cubed(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(1,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    i = randint(0,number_of_particles-2)
    j = randint(i+1,number_of_particles-1)
    r_ij = make_r(x)(i,j)
    if number_of_dimensions == 1:
        k = 0
    else:
        k = randint(0,number_of_dimensions-1)

    substitutions = generate_random_substitutions(*x.ravel())
    self.assertAlmostEqual((r_ij**3).diff(x[i,k]).subs(substitutions),(3.0*r_ij*(+(x[i,k]-x[j,k]))).subs(substitutions))
    self.assertAlmostEqual((r_ij**3).diff(x[j,k]).subs(substitutions),(3.0*r_ij*(-(x[i,k]-x[j,k]))).subs(substitutions))
#@-node:gcross.20090906145419.1345:test_diff_r_cubed
#@+node:gcross.20090904201537.1563:test_diff_rho
@with_checker
def test_diff_rho(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    i = randint(0,number_of_particles-1)
    rho_i = make_rho(x)(i)

    substitutions = generate_random_substitutions(*x.ravel())

    for k in xrange(number_of_dimensions):
        if k == 0 or k == 1:
            self.assertAlmostEqual(rho_i.diff(x[i,k]).subs(substitutions),(x[i,k]/rho_i).subs(substitutions))
            ip = i
            while ip == i:
                ip = randint(0,number_of_particles-1)
            self.assertAlmostEqual(rho_i.diff(x[ip,k]).subs(substitutions),0)
        else:
            for ip in xrange(number_of_particles):
                self.assertAlmostEqual(rho_i.diff(x[ip,k]).subs(substitutions),0)
#@-node:gcross.20090904201537.1563:test_diff_rho
#@+node:gcross.20090906145419.1351:test_diff_rho_squared
@with_checker
def test_diff_rho_squared(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    i = randint(0,number_of_particles-1)
    rho_i = make_rho(x)(i)

    substitutions = generate_random_substitutions(*x.ravel())

    rho_i_squared = rho_i**2

    for k in xrange(number_of_dimensions):
        if k == 0 or k == 1:
            self.assertAlmostEqual(rho_i_squared.diff(x[i,k]).subs(substitutions),(2.0*x[i,k]).subs(substitutions))
            ip = i
            while ip == i:
                ip = randint(0,number_of_particles-1)
            self.assertAlmostEqual(rho_i_squared.diff(x[ip,k]).subs(substitutions),0)
        else:
            for ip in xrange(number_of_particles):
                self.assertAlmostEqual(rho_i_squared.diff(x[ip,k]).subs(substitutions),0)
#@-node:gcross.20090906145419.1351:test_diff_rho_squared
#@+node:gcross.20090907141528.1353:test_gradient_backflow_numerator
@with_checker(number_of_calls=100)
def test_gradient_backflow_numerator(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    for i in xrange(number_of_particles):
        for j in xrange(number_of_particles):
            if i == j: continue
            numerator = backflow_numerator(x,i,j)
            for ip in xrange(number_of_particles):
                for kp in xrange(number_of_dimensions):
                    computed = numerator.diff(x[ip,kp])
                    formula = gradient_backflow_numerator(a,r,rho,x,i,j,ip,kp)
                    self.assertAlmostEqual(computed.subs(substitutions),formula.subs(substitutions))
#@-node:gcross.20090907141528.1353:test_gradient_backflow_numerator
#@+node:gcross.20090905201409.1340:test_gradient_backflow_denominator
@with_checker(number_of_calls=100)
def test_gradient_backflow_denominator(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    substitutions = generate_random_substitutions(a,*(x.ravel()))


    for i in xrange(number_of_particles):
        for j in xrange(number_of_particles):
            if i == j: continue
            denominator = backflow_denominator(a,r,rho,x,i,j)
            for ip in xrange(number_of_particles):
                for kp in xrange(number_of_dimensions):
                    computed = denominator.diff(x[ip,kp])
                    formula = gradient_backflow_denominator(a,r,rho,x,i,j,ip,kp)
                    self.assertAlmostEqual(computed.subs(substitutions),formula.subs(substitutions))
#@-node:gcross.20090905201409.1340:test_gradient_backflow_denominator
#@+node:gcross.20090907141528.1357:test_gradient_backflow_via_ratio_rule
@with_checker(number_of_calls=100)
def test_gradient_backflow_via_ratio_rule(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    b = backflow(a,r,rho,x)

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    for ip in xrange(number_of_particles):
        for kp in xrange(number_of_dimensions):
            computed = b.diff(x[ip,kp])

            formula = 0

            for i in xrange(number_of_particles):
                for j in xrange(number_of_particles):
                    if i == j: continue
                    numerator = backflow_numerator(x,i,j)
                    denominator = backflow_denominator(a,r,rho,x,i,j)
                    grad_numerator = gradient_backflow_numerator(a,r,rho,x,i,j,ip,kp)
                    grad_denominator = gradient_backflow_denominator(a,r,rho,x,i,j,ip,kp)

                    formula += grad_numerator/denominator - numerator/denominator**2 * grad_denominator

            formula *= a

            self.assertAlmostEqual(computed.subs(substitutions),formula.subs(substitutions))
#@-node:gcross.20090907141528.1357:test_gradient_backflow_via_ratio_rule
#@+node:gcross.20090907141528.1365:test_gradient_backflow_numerator_resummation
@with_checker(number_of_calls=100)
def test_gradient_backflow_numerator_resummation(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    b = backflow(a,r,rho,x)

    def D(i,j):
        return r(i,j)**3 * (1-a/r(i,j))

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    for ip in xrange(number_of_particles):
        for kp in xrange(number_of_dimensions):
            sum_over_ij = 0
            for i in xrange(number_of_particles):
                for j in xrange(number_of_particles):
                    if i == j: continue
                    grad_numerator = gradient_backflow_numerator(a,r,rho,x,i,j,ip,kp)
                    denominator = backflow_denominator(a,r,rho,x,i,j)

                    sum_over_ij += grad_numerator/denominator

            rho_ip = rho(ip)

            d_kp0 = 1 if kp == 0 else 0
            d_kp1 = 1 if kp == 1 else 0

            sum_over_j_only = 0
            for j in xrange(number_of_particles):
                if ip == j: continue
                sum_over_j_only += 1.0/D(ip,j)*(x[j,1]*d_kp0-x[j,0]*d_kp1)*(1/rho_ip**2 - 1/rho(j)**2)

            self.assertAlmostEqual(sum_over_ij.subs(substitutions),sum_over_j_only.subs(substitutions))
#@-node:gcross.20090907141528.1365:test_gradient_backflow_numerator_resummation
#@+node:gcross.20090907141528.1371:test_gradient_backflow_denominator_over_denominator
@with_checker(number_of_calls=100)
def test_gradient_backflow_denominator_over_denominator(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    substitutions = generate_random_substitutions(a,*(x.ravel()))


    for i in xrange(number_of_particles):
        for j in xrange(number_of_particles):
            if i == j: continue
            denominator = backflow_denominator(a,r,rho,x,i,j)

            rho_i = rho(i)
            r_ij = r(i,j)

            for ip in xrange(number_of_particles):
                d_ipi = 1 if ip == i else 0
                d_ipj = 1 if ip == j else 0
                for kp in xrange(number_of_dimensions):
                    d_kp0 = 1 if kp == 0 else 0
                    d_kp1 = 1 if kp == 1 else 0
                    computed = gradient_backflow_denominator(a,r,rho,x,i,j,ip,kp)/denominator
                    formula = 2*x[ip,kp]*d_ipi*(d_kp0+d_kp1)/rho_i**2 + 3*(x[i,kp]-x[j,kp])*(d_ipi-d_ipj)/r_ij**2 + a*(x[i,kp]-x[j,kp])*(d_ipi-d_ipj)/(r_ij**3*(1-a/r_ij))
                    self.assertAlmostEqual(computed.subs(substitutions),formula.subs(substitutions))
#@-node:gcross.20090907141528.1371:test_gradient_backflow_denominator_over_denominator
#@+node:gcross.20090907141528.1367:test_gradient_backflow_denominator_rewrite_1
@with_checker(number_of_calls=100)
def test_gradient_backflow_denominator_rewrite_1(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    b = backflow(a,r,rho,x)

    def C(i,j):
        return x[i,0]*x[j,1]-x[i,1]*x[j,0]

    def D(i,j):
        return r(i,j)**3 * (1-a/r(i,j))

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    for ip in xrange(number_of_particles):
        for kp in xrange(number_of_dimensions):
            sum_1 = 0
            for i in xrange(number_of_particles):
                for j in xrange(number_of_particles):
                    if i == j: continue
                    numerator = backflow_numerator(x,i,j)
                    denominator = backflow_denominator(a,r,rho,x,i,j)
                    grad_denominator = gradient_backflow_denominator(a,r,rho,x,i,j,ip,kp)

                    sum_1 += numerator/denominator**2 * grad_denominator

            rho_ip = rho(ip)

            d_kp0 = 1 if kp == 0 else 0
            d_kp1 = 1 if kp == 1 else 0

            sum_2 = 0
            for i in xrange(number_of_particles):
                for j in xrange(number_of_particles):
                    d_ipi = 1 if ip == i else 0
                    d_ipj = 1 if ip == j else 0
                    if i == j: continue
                    r_ij = r(i,j)
                    rho_i = rho(i)
                    sum_2 += (C(i,j)/(rho_i**2*D(i,j))) * ((2*x[ip,kp]*d_ipi*(d_kp0+d_kp1))/rho_i**2 + 3*(x[i,kp]-x[j,kp])*(d_ipi-d_ipj)/r_ij**2 + a*(x[i,kp]-x[j,kp])*(d_ipi-d_ipj)/(r_ij**3*(1-a/r_ij)))

            self.assertAlmostEqual(sum_1.subs(substitutions),sum_2.subs(substitutions))
#@-node:gcross.20090907141528.1367:test_gradient_backflow_denominator_rewrite_1
#@+node:gcross.20090907141528.1373:test_gradient_backflow_denominator_rewrite_2
@with_checker(number_of_calls=100)
def test_gradient_backflow_denominator_rewrite_2(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    b = backflow(a,r,rho,x)

    def C(i,j):
        return x[i,0]*x[j,1]-x[i,1]*x[j,0]

    def D(i,j):
        return r(i,j)**3 * (1-a/r(i,j))

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    for ip in xrange(number_of_particles):
        for kp in xrange(number_of_dimensions):
            sum_1 = 0
            for i in xrange(number_of_particles):
                for j in xrange(number_of_particles):
                    if i == j: continue
                    numerator = backflow_numerator(x,i,j)
                    denominator = backflow_denominator(a,r,rho,x,i,j)
                    grad_denominator = gradient_backflow_denominator(a,r,rho,x,i,j,ip,kp)

                    sum_1 += numerator/denominator**2 * grad_denominator

            rho_ip = rho(ip)

            d_kp0 = 1 if kp == 0 else 0
            d_kp1 = 1 if kp == 1 else 0

            sum_2 = 0
            for i in xrange(number_of_particles):
                for j in xrange(number_of_particles):
                    d_ipi = 1 if ip == i else 0
                    d_ipj = 1 if ip == j else 0
                    if i == j: continue
                    r_ij = r(i,j)
                    rho_i = rho(i)
                    C_ij = C(i,j)
                    D_ij = D(i,j)
                    sum_2 += 2*C_ij*x[ip,kp]*d_ipi*(d_kp0+d_kp1)/D_ij/rho_i**4 + 3*C_ij*(x[i,kp]-x[j,kp])*(d_ipi-d_ipj)/D_ij/rho_i**2/r_ij**2 + C_ij*a*(x[i,kp]-x[j,kp])*(d_ipi-d_ipj)/D_ij/rho_i**2/r_ij**3/(1-a/r_ij)

            self.assertAlmostEqual(sum_1.subs(substitutions),sum_2.subs(substitutions))
#@-node:gcross.20090907141528.1373:test_gradient_backflow_denominator_rewrite_2
#@+node:gcross.20090907214729.1357:test_gradient_backflow_denominator_rewrite_3
@with_checker(number_of_calls=100)
def test_gradient_backflow_denominator_rewrite_3(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    b = backflow(a,r,rho,x)

    def C(i,j):
        return x[i,0]*x[j,1]-x[i,1]*x[j,0]

    def D(i,j):
        return r(i,j)**3 * (1-a/r(i,j))

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    for ip in xrange(number_of_particles):
        for kp in xrange(number_of_dimensions):
            sum_1 = 0
            for i in xrange(number_of_particles):
                for j in xrange(number_of_particles):
                    if i == j: continue
                    numerator = backflow_numerator(x,i,j)
                    denominator = backflow_denominator(a,r,rho,x,i,j)
                    grad_denominator = gradient_backflow_denominator(a,r,rho,x,i,j,ip,kp)

                    sum_1 += numerator/denominator**2 * grad_denominator

            rho_ip = rho(ip)

            d_kp0 = 1 if kp == 0 else 0
            d_kp1 = 1 if kp == 1 else 0

            sum_2 = 0
            for j in xrange(number_of_particles):
                if j == ip: continue
                sum_2 += (
                    2*C(ip,j)*x[ip,kp]*(d_kp0+d_kp1)/D(ip,j)/rho(ip)**4
                  + 3*C(ip,j)*(x[ip,kp]-x[j,kp])/D(ip,j)/rho(ip)**2/r(ip,j)**2 - 3*C(j,ip)*(x[j,kp]-x[ip,kp])/D(j,ip)/rho(j)**2/r(j,ip)**2
                  + C(ip,j)*a*(x[ip,kp]-x[j,kp])/D(ip,j)/rho(ip)**2/r(ip,j)**3/(1-a/r(ip,j)) - C(j,ip)*a*(x[j,kp]-x[ip,kp])/D(j,ip)/rho(j)**2/r(j,ip)**3/(1-a/r(j,ip))
                )

            self.assertAlmostEqual(sum_1.subs(substitutions),sum_2.subs(substitutions))
#@-node:gcross.20090907214729.1357:test_gradient_backflow_denominator_rewrite_3
#@+node:gcross.20090907141528.1369:test_gradient_backflow_denominator_resummation
@with_checker(number_of_calls=100)
def test_gradient_backflow_denominator_resummation(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    b = backflow(a,r,rho,x)

    def C(i,j):
        return x[i,0]*x[j,1]-x[i,1]*x[j,0]

    def D(i,j):
        return r(i,j)**3 * (1-a/r(i,j))

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    for ip in xrange(number_of_particles):
        for kp in xrange(number_of_dimensions):
            sum_over_ij = 0
            for i in xrange(number_of_particles):
                for j in xrange(number_of_particles):
                    if i == j: continue
                    numerator = backflow_numerator(x,i,j)
                    denominator = backflow_denominator(a,r,rho,x,i,j)
                    grad_denominator = gradient_backflow_denominator(a,r,rho,x,i,j,ip,kp)

                    sum_over_ij += numerator/denominator**2 * grad_denominator

            rho_ip = rho(ip)

            d_kp0 = 1 if kp == 0 else 0
            d_kp1 = 1 if kp == 1 else 0

            sum_over_j_only = 0
            for j in xrange(number_of_particles):
                if ip == j: continue
                r_ipj = r(ip,j)
                sum_over_j_only += (C(ip,j)/D(ip,j))*(2*x[ip,kp]*(d_kp0 + d_kp1)/rho_ip**4 + (x[ip,kp]-x[j,kp])/r_ipj**2 * (1/rho_ip**2-1/rho(j)**2) * (3+a/(r_ipj*(1-a/r_ipj))))

            self.assertAlmostEqual(sum_over_ij.subs(substitutions),sum_over_j_only.subs(substitutions))
#@-node:gcross.20090907141528.1369:test_gradient_backflow_denominator_resummation
#@+node:gcross.20090907141528.1345:test_gradient_backflow
@with_checker(number_of_calls=10)
def test_gradient_backflow(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
  ):
    x = make_x(number_of_particles,number_of_dimensions)
    a = Symbol('a')
    r = make_r(x)
    rho = make_rho(x)

    b = backflow(a,r,rho,x)

    substitutions = generate_random_substitutions(a,*(x.ravel()))

    for ip in xrange(number_of_particles):
        for kp in xrange(number_of_dimensions):
            computed = b.diff(x[ip,kp])
            formula = backflow_gradient(a,r,rho,x,ip,kp)
            self.assertAlmostEqual(computed.subs(substitutions),formula.subs(substitutions))
#@-node:gcross.20090907141528.1345:test_gradient_backflow
#@+node:gcross.20090908092630.1348:test_fortran_gradient_backflow
@with_checker(number_of_calls=100)
def test_fortran_gradient_backflow(self,
    number_of_particles=irange(2,5),
    number_of_dimensions=irange(2,5),
    hard_sphere_radius=unit_interval_float,
    rotation_rate=unit_interval_float,
  ):
    def make_compute_formula_gradient_backflow(a_value,x_values):
        x = make_x(number_of_particles,number_of_dimensions)
        a = Symbol('a')
        r = make_r(x)
        rho = make_rho(x)
        substitutions = dict(itertools.izip(x.ravel(),x_values.ravel()))
        substitutions[a] = a_value

        def compute_formula_gradient_backflow(ip,kp):
            formula = backflow_gradient(a,r,rho,x,ip,kp)
            return formula.subs(substitutions)

        return compute_formula_gradient_backflow

    x = rand(number_of_particles,number_of_dimensions)
    compute_formula_gradient_backflow = make_compute_formula_gradient_backflow(hard_sphere_radius,x)

    x = array(x.reshape((1,)+x.shape),dtype=double,order='Fortran')
    xij2 = zeros((1,number_of_particles,number_of_particles),dtype=double,order='Fortran')
    vpif.xij.update_xij(xij2,x)

    gradient_backflow = vpif.hard_sphere_interaction.compute_gradient_backflow(x,xij2,hard_sphere_radius,rotation_rate,1,2)

    for ip in xrange(number_of_particles):
        for kp in xrange(number_of_dimensions):
            computed = gradient_backflow[0,ip,kp]
            formula = compute_formula_gradient_backflow(ip,kp)
            self.assertAlmostEqual(computed,rotation_rate*formula)
#@-node:gcross.20090908092630.1348:test_fortran_gradient_backflow
#@-node:gcross.20090904201537.1795:(verified)
#@+node:gcross.20090907141528.1361:(unverified)
#@-node:gcross.20090907141528.1361:(unverified)
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
