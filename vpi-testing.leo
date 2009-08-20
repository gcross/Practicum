<?xml version="1.0" encoding="utf-8"?>
<?xml-stylesheet ekr_test?>
<leo_file>
<leo_header file_format="2" tnodes="0" max_tnode_index="0" clone_windows="0"/>
<globals body_outline_ratio="0.5">
	<global_window_position top="10" left="432" height="731" width="738"/>
	<global_log_window_position top="0" left="0" height="0" width="0"/>
</globals>
<preferences/>
<find_panel_settings/>
<vnodes>
<v t="gcross.20090817094245.1729" str_leo_pos="7,2"><vh>@chapters</vh></v>
<v t="gcross.20090817100500.1247"><vh>@button Close plot</vh></v>
<v t="gcross.20090818083714.3694" a="E"><vh>@settings</vh>
<v t="gcross.20090818083714.3695"><vh>@string ipython_argv = ipython -pylab</vh></v>
</v>
<v t="gcross.20090817094840.1236" a="E"><vh>@ipy-startup</vh>
<v t="gcross.20090817094840.1237"><vh>import vpi</vh></v>
<v t="gcross.20090817094840.1241"><vh>import numpy</vh></v>
<v t="gcross.20090817100500.1362"><vh>import pylab</vh></v>
<v t="gcross.20090820084503.1342"><vh>import scipy</vh></v>
<v t="gcross.20090820084503.1343"><vh>import itertools</vh></v>
<v t="gcross.20090820084503.1428"><vh>import functools</vh></v>
<v t="gcross.20090820145419.1287"><vh>import math.atan</vh></v>
<v t="gcross.20090817094840.1240"><vh>make_string</vh></v>
</v>
<v t="gcross.20090817094840.1235"><vh>Console</vh></v>
<v t="gcross.20090820084503.1255"><vh>Project: test U_eff</vh>
<v t="gcross.20090817094840.1239"><vh>compute_effective_rotational_potential documentation</vh></v>
<v t="gcross.20090817102809.1253"><vh>computer_angular_derivatives</vh></v>
<v t="gcross.20090817100924.1251"><vh>Scan derivatives</vh></v>
<v t="gcross.20090817104911.1258"><vh>Scan clumping</vh></v>
<v t="gcross.20090817094840.1238"><vh>Scan potentials</vh></v>
</v>
<v t="gcross.20090820084503.1256" a="E"><vh>Project: test derivatives</vh>
<v t="gcross.20090820084503.1257"><vh>phase</vh></v>
<v t="gcross.20090820084503.1259"><vh>make_phase1</vh></v>
<v t="gcross.20090820084503.1426"><vh>make_phase2</vh></v>
<v t="gcross.20090820084503.1258"><vh>data</vh></v>
<v t="gcross.20090820084503.1341"><vh>first derivatives</vh></v>
<v t="gcross.20090820084503.1427"><vh>second derivatives</vh></v>
<v t="gcross.20090820145419.1288"><vh>my code</vh></v>
</v>
<v t="gcross.20090820084503.1429" a="E"><vh>Project: test grad U_eff</vh>
<v t="gcross.20090820084503.1430"><vh>def U</vh></v>
<v t="gcross.20090820084503.1433"><vh>make_U</vh></v>
<v t="gcross.20090820084503.1432"><vh>data</vh></v>
<v t="gcross.20090820084503.1435"><vh>derivatives</vh></v>
</v>
</vnodes>
<tnodes>
<t tx="gcross.20090817094245.1729"></t>
<t tx="gcross.20090817094840.1235">print potentials
</t>
<t tx="gcross.20090817094840.1236"></t>
<t tx="gcross.20090817094840.1237">import vpi
</t>
<t tx="gcross.20090817094840.1238">del vpi
import vpi
print type(vpi.angular_momentum)
compute_effective_rotational_potential = vpi.angular_momentum.compute_effective_rotational_potential

n_particles = 10
n_slices = 1
n_dimensions = 3

x = array(rand(n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')
for i in xrange(n_slices):
    for j in xrange(n_particles):
        x[i,j,2] = 0;
        x[i,j,:] /= norm(x[i,j,:])

plane_axis_1 = 1
plane_axis_2 = 2

move_start = 1
move_end = 1

U = zeros((n_slices,n_particles),dtype=double,order='Fortran')
gradU = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')

close()
for frame_angular_velocity in arange(0,1,0.1):
    potentials = []
    for fixed_angular_momentum in xrange(2,n_particles):
        U[...] = 0
        compute_effective_rotational_potential(x,plane_axis_1,plane_axis_2,frame_angular_velocity,fixed_angular_momentum,move_start,move_end,U,gradU)
        potentials.append((fixed_angular_momentum,sum(U[0])))
    X, Y = zip(*potentials)
    print frame_angular_velocity, argmin(Y)+2
    plot(X,Y)
</t>
<t tx="gcross.20090817094840.1239">@cl make_string


compute_effective_rotational_potential - Function signature:
  compute_effective_rotational_potential(x,fixed_rotation_axis,frame_angular_velocity,fixed_angular_momentum,move_start,move_end,U_rot,[n_slices,n_particles,n_dimensions])
Required arguments:
  x : input rank-3 array('d') with bounds (n_slices,n_particles,n_dimensions)
  fixed_rotation_axis : input int
  frame_angular_velocity : input float
  fixed_angular_momentum : input int
  move_start : input int
  move_end : input int
  U_rot : in/output rank-2 array('d') with bounds (n_slices,n_particles)
Optional arguments:
  n_slices := shape(x,0) input int
  n_particles := shape(x,1) input int
  n_dimensions := shape(x,2) input int
</t>
<t tx="gcross.20090817094840.1240">def make_string(body,_):
    return body
</t>
<t tx="gcross.20090817094840.1241">from numpy import *
from numpy.linalg import *
</t>
<t tx="gcross.20090817100500.1247">import pylab
pylab.close()
</t>
<t tx="gcross.20090817100500.1362">from pylab import *
</t>
<t tx="gcross.20090817100924.1251">del vpi
import vpi
print type(vpi.angular_momentum)
compute_angular_derivatives = vpi.angular_momentum.compute_angular_derivatives

n_particles = 100
n_dimensions = 3

x = array(rand(n_particles,n_dimensions),dtype=double,order='Fortran')

fixed_rotation_axis = 3

close()
potentials = []
for fixed_angular_momentum in xrange(0,n_particles):
    derivatives = compute_angular_derivatives(x,fixed_rotation_axis,fixed_angular_momentum)
    potentials.append((fixed_angular_momentum,sum(derivatives*(derivatives-2))))
X, Y = zip(*potentials)
plot(X,Y)

</t>
<t tx="gcross.20090817102809.1253">@cl make_string


compute_angular_derivatives - Function signature:
  derivatives = compute_angular_derivatives(x,fixed_rotation_axis,fixed_angular_momentum,[N_particles,N_dimensions])
Required arguments:
  x : input rank-2 array('d') with bounds (N_particles,N_dimensions)
  fixed_rotation_axis : input int
  fixed_angular_momentum : input int
Optional arguments:
  N_particles := shape(x,0) input int
  N_dimensions := shape(x,1) input int
Return objects:
  derivatives : rank-1 array('d') with bounds (N_particles)</t>
<t tx="gcross.20090817104911.1258">del vpi
import vpi
print type(vpi.angular_momentum)
compute_effective_rotational_potential = vpi.angular_momentum.compute_effective_rotational_potential

n_particles = 10
n_slices = 1
n_dimensions = 3

x = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')

fixed_rotation_axis = 3

frame_angular_velocity = 0

move_start = 1
move_end = 1

U_rot = zeros((n_slices,n_particles),dtype=double,order='Fortran')

close()
for fixed_angular_momentum in xrange(1,n_particles):
    potentials = []
    for angular_width in arange(0,0.7,0.01):
        angles = array(range(n_particles))*angular_width
        radii = rand(n_slices,n_particles)
        x[:,:,0] = radii*cos(angles)
        x[:,:,1] = radii*sin(angles)
        x[:,:,2] = rand(n_slices,n_particles)
        U_rot[...] = 0
        compute_effective_rotational_potential(x,fixed_rotation_axis,frame_angular_velocity,fixed_angular_momentum,move_start,move_end,U_rot)
        potentials.append((angular_width,sum(U_rot[0])))
    X, Y = zip(*potentials)
    print X, Y
    plot(X,Y)
</t>
<t tx="gcross.20090818083714.3694"></t>
<t tx="gcross.20090818083714.3695"></t>
<t tx="gcross.20090820084503.1255"></t>
<t tx="gcross.20090820084503.1256">@others
</t>
<t tx="gcross.20090820084503.1257">def phase(N_rotating_particles, angles):
    C = 0
    S = 0
    for a in combinations(angles,N_rotating_particles):
        C += cos(sum(a))
        S += sin(sum(a))
    return atan(S/C)
</t>
<t tx="gcross.20090820084503.1258">N_particles = 4
N_rotating_particles = 1
N_dimensions = 5

x = rand(N_particles,N_dimensions)
angles = arctan2(x[:,1],x[:,0])
</t>
<t tx="gcross.20090820084503.1259">def make_phase1(N_rotating_particles,angles,index):
    angles = angles.copy()
    def phase1(angle):
        angles[index] = angle
        return phase(N_rotating_particles,angles)
    return phase1
</t>
<t tx="gcross.20090820084503.1341">diffs = [derivative(make_phase1(N_rotating_particles,angles,i),angles[i],dx=1e-6,n=1,order=13) for i in xrange(len(angles))]
print diffs
</t>
<t tx="gcross.20090820084503.1342">from scipy.misc import derivative
</t>
<t tx="gcross.20090820084503.1343">from itertools import *</t>
<t tx="gcross.20090820084503.1426">def make_phase2(N_rotating_particles,angles,index1,index2):
    angles = angles.copy()
    def phase2(angle1,angle2):
        angles[index1] = angle1
        angles[index2] = angle2
        return phase(N_rotating_particles,angles)
    return phase2
</t>
<t tx="gcross.20090820084503.1427">def phase_d2(N_rotating_particles,angles,index1,index2):
    phase2 = make_phase2(2,angles,index1,index2)
    def diff_phase(angle1):
        phase1 = partial(phase2,angle1)
        return derivative(phase1,angles[index2],dx=1e-6,n=1,order=13)
    return derivative(diff_phase,angles[index1],dx=1e-6,n=1,order=13)

diff2 = zeros((len(angles),len(angles)),dtype=double)
for i in xrange(len(angles)):
    for j in xrange(len(angles)):
        if i == j:
            diff2[i,i] = derivative(make_phase1(2,angles,i),angles[i],dx=1e-6,n=2,order=13)
        else:
            diff2[i,j] = phase_d2(2,angles,i,j)
print diff2

</t>
<t tx="gcross.20090820084503.1428">from functools import *
</t>
<t tx="gcross.20090820084503.1429">@others
</t>
<t tx="gcross.20090820084503.1430">def compute_Ueff(x):
    U = zeros(x.shape[:2],dtype=double,order='Fortran')
    gradU = zeros(x.shape,dtype=double,order='Fortran')
    vpi.angular_momentum.compute_effective_rotational_potential(
        x, lambda_,
        1, 2, frame_angular_velocity, N_rotating_particles,
        U, gradU
    )
    return sum(U)

def compute_gradUeff(x):
    U = zeros(x.shape[:2],dtype=double,order='Fortran')
    gradU = zeros(x.shape,dtype=double,order='Fortran')
    vpi.angular_momentum.compute_effective_rotational_potential(
        x, lambda_,
        1, 2, frame_angular_velocity, N_rotating_particles,
        U, gradU
    )
    return gradU
</t>
<t tx="gcross.20090820084503.1432">N_slices = 1
N_particles = 9
N_dimensions = 3

lambda_ = 1.3
frame_angular_velocity = 0.4
N_rotating_particles = 1

x = rand(N_slices,N_particles,N_dimensions)*10-5
</t>
<t tx="gcross.20090820084503.1433">def make_U(x,slice_,particle,coordinate):
    x = x.copy()
    def U(v):
        x[slice_,particle,coordinate] = v
        return compute_Ueff(x)
    return U
</t>
<t tx="gcross.20090820084503.1435">diffs = array(
    [[[
        derivative(make_U(x,i,j,k),x[i,j,k],dx=1e-6,n=1,order=13)
        for k in xrange(x.shape[-1])]
        for j in xrange(x.shape[-2])]
        for i in xrange(x.shape[-3])]
    )
print diffs
print compute_gradUeff(x)
</t>
<t tx="gcross.20090820145419.1287">from math import atan
</t>
<t tx="gcross.20090820145419.1288">vpi.angular_momentum.compute_angular_derivatives(x,1,2,N_rotating_particles)
</t>
</tnodes>
</leo_file>
