#@+leo-ver=4-thin
#@+node:gcross.20090819152718.1279:@thin run.py
import sys
sys.path.append("../../lib")

#@<< Imports >>
#@+node:gcross.20090819152718.1280:<< Imports >>
import vpi

from numpy import *
from numpy.random import rand

from mpi4py import MPI

comm = MPI.COMM_WORLD

number_of_processors = comm.Get_size()
my_rank = comm.Get_rank()
#@nonl
#@-node:gcross.20090819152718.1280:<< Imports >>
#@nl

#@<< Configuration >>
#@+node:gcross.20090819152718.1281:<< Configuration >>
# System parameters
n_particles = 2
n_slices = 102
lambda_ = 0.5
n_dimensions = 3

# Angular momentum parameters
frame_angular_velocity = 1
N_rotating_particles = 1
rotation_plane_axis_1 = 1
rotation_plane_axis_2 = 2

# Run parameters
total_number_of_observations = 100000
number_of_prethermalization_steps = 1000

# Histogram parameters
_1d_densities_histogram_left  = [-2,-2,-2]
_1d_densities_histogram_right = [+2,+2,+2]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

# Move parameters
dM = 22
move_type_probabilities = [0.9,0.1,0]
move_type_differentials = [0.1,0.6,0]
low_swap_dimension = 1
high_swap_dimension= 3

# Potential parameters
harmonic_oscillator_coefficients = array((1,1,1),dtype=double)

# Miscellaneous
use_4th_order_green_function = True
#@-node:gcross.20090819152718.1281:<< Configuration >>
#@nl

#@+others
#@+node:gcross.20090819152718.1282:Functions
#@+node:gcross.20090819152718.1283:Physics
#@+node:gcross.20090819152718.1284:Potential function
def compute_potential(x,xij2):
    x_sq = x**2
    U = array(dot(x_sq,harmonic_oscillator_coefficients)/2.0,dtype=double,order='Fortran')
    gradU  = array(x,dtype=double,order='Fortran')
    gradU *= harmonic_oscillator_coefficients
    vpi.angular_momentum.compute_effective_rotational_potential(
        x,lambda_,
        rotation_plane_axis_1,rotation_plane_axis_2,
        frame_angular_velocity,N_rotating_particles,
        U, gradU
    )
    gradU2 = sum(sum(gradU**2,axis=-1),axis=-1)
    return U, gradU2, False
#@-node:gcross.20090819152718.1284:Potential function
#@+node:gcross.20090819152718.1285:Trial function
def trial_function(x,_):
    return -sum(dot(x**2,harmonic_oscillator_coefficients))/2
#@-node:gcross.20090819152718.1285:Trial function
#@+node:gcross.20090819152718.1286:Trial derivatives
def trial_derivatives(x,xij2):
    gradient_of_log_trial_fn = -x*harmonic_oscillator_coefficients
    laplacian_of_log_trial_fn = -sum(harmonic_oscillator_coefficients)*x.shape[0]
    return gradient_of_log_trial_fn, laplacian_of_log_trial_fn
#@-node:gcross.20090819152718.1286:Trial derivatives
#@-node:gcross.20090819152718.1283:Physics
#@+node:gcross.20090819152718.1287:Observables
#@+node:gcross.20090819152718.1288:Energy
def compute_energy_estimates():
    estimates = []
    for i in [0,-1]:
        gradient_of_log_trial_fn, laplacian_of_log_trial_fn = trial_derivatives(x[i],xij2[i])
        estimates.append(
            vpi.observables.compute_local_energy_estimate(
                U[i],
                gradient_of_log_trial_fn, laplacian_of_log_trial_fn,
                lambda_,
            )
        )
    return estimates
#@-node:gcross.20090819152718.1288:Energy
#@-node:gcross.20090819152718.1287:Observables
#@+node:gcross.20090819152718.1289:Histograms
#@+node:gcross.20090819152718.1290:write_histogram
def write_histogram(filename,histogram,left,right):
    with open(filename,"w") as f:
        bin_count = len(histogram)
        bin_width = float(right-left)/_1d_densities_histogram_bin_count
        current = float(left)
        for count in histogram:
            print >> f, "{0} {1}".format(current+bin_width/2,float(count)/(total_number_of_observations*n_particles))
            current += bin_width
#@-node:gcross.20090819152718.1290:write_histogram
#@+node:gcross.20090819152718.1291:write_radial_histogram
def write_radial_histogram(filename,histogram,radius):
    with open(filename,"w") as f:
        bin_count = len(histogram)
        bin_width = float(radius)/radial_densities_histogram_bin_count
        current = 0
        for count in histogram:
            normalization = 4.0/3.0*pi*((current+bin_width)**3-current**3)
            print >> f, "{0} {1}".format(current+bin_width/2,float(count)/(total_number_of_observations*n_particles*normalization))     
            current += bin_width
#@-node:gcross.20090819152718.1291:write_radial_histogram
#@+node:gcross.20090819152718.1292:compute_total_histogram
def compute_total_histogram(histogram):
    total_histogram = zeros(histogram.shape,dtype='i',order='Fortran')
    comm.Reduce((histogram,MPI.INT),(total_histogram,MPI.INT))
    return total_histogram
#@-node:gcross.20090819152718.1292:compute_total_histogram
#@-node:gcross.20090819152718.1289:Histograms
#@-node:gcross.20090819152718.1282:Functions
#@-others

#@<< Initialization >>
#@+node:gcross.20090819152718.1293:<< Initialization >>
vpi.rand_utils.init_seed(my_rank)

n_thermalizations_per_observation = n_particles*n_slices//dM

assert (n_slices % 2 == 0 and n_slices % 4 == 2)
c_slice = n_slices / 2

x = vpi.lattice.make_lattice(1.0,n_slices,n_particles,n_dimensions)

xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
vpi.xij.update_xij(xij2,x)

U = zeros((n_slices,n_particles),dtype=double,order='Fortran')
gradU2 = zeros((n_slices),dtype=double,order='Fortran')

slice_move_attempted_counts = zeros((n_slices,),'i')
slice_move_accepted_counts = zeros((n_slices,),'i')

move_type_attempted_counts = zeros((3,),'i')
move_type_accepted_counts = zeros((3,),'i')

U_weights, gU2_weights = vpi.gfn.initialize_4th_order_weights(n_slices)

_1d_densities_histogram_left = array(_1d_densities_histogram_left,dtype=double,order='Fortran')
_1d_densities_histogram_right = array(_1d_densities_histogram_right,dtype=double,order='Fortran')

left_slice_1d_densities_histogram = zeros((n_dimensions,_1d_densities_histogram_bin_count),dtype='i',order='Fortran')
center_slice_1d_densities_histogram = zeros((n_dimensions,_1d_densities_histogram_bin_count),dtype='i',order='Fortran')
right_slice_1d_densities_histogram = zeros((n_dimensions,_1d_densities_histogram_bin_count),dtype='i',order='Fortran')

left_slice_radial_densities_histogram = zeros((radial_densities_histogram_bin_count),dtype='i',order='Fortran')
center_slice_radial_densities_histogram = zeros((radial_densities_histogram_bin_count),dtype='i',order='Fortran')
right_slice_radial_densities_histogram = zeros((radial_densities_histogram_bin_count),dtype='i',order='Fortran')
#@-node:gcross.20090819152718.1293:<< Initialization >>
#@nl

#@<< Main iteration >>
#@+node:gcross.20090819152718.1294:<< Main iteration >>
vpi.thermalize.thermalize_path(
    x,xij2,
    U,gradU2,
    number_of_prethermalization_steps,
    move_type_probabilities,move_type_differentials,
    dM,
    lambda_,
    low_swap_dimension,high_swap_dimension,
    slice_move_attempted_counts,move_type_attempted_counts,
    slice_move_accepted_counts,move_type_accepted_counts,
    compute_potential,trial_function,
    U_weights,gU2_weights,
    use_4th_order_green_function,
)

energy_estimates = zeros((2,),dtype=double)

number_of_observations = total_number_of_observations // number_of_processors + 1

total_number_of_observations = number_of_observations * number_of_processors

for _ in xrange(number_of_observations):
    vpi.thermalize.thermalize_path(
        x,xij2,
        U,gradU2,
        n_thermalizations_per_observation,
        move_type_probabilities,move_type_differentials,
        dM,
        lambda_,
        low_swap_dimension,high_swap_dimension,
        slice_move_attempted_counts,move_type_attempted_counts,
        slice_move_accepted_counts,move_type_accepted_counts,
        compute_potential,trial_function,
        U_weights,gU2_weights,
        use_4th_order_green_function,
    )
    energy_estimates += compute_energy_estimates()
    for slice_number, _1d_densities_histogram, radial_densities_histogram in [
            (0,left_slice_1d_densities_histogram,left_slice_radial_densities_histogram),
            (c_slice-1,center_slice_1d_densities_histogram,center_slice_radial_densities_histogram),
            (-1,right_slice_1d_densities_histogram,right_slice_radial_densities_histogram),
        ]:
        vpi.histograms.accumulate_1d_densities(
            x[slice_number],
            _1d_densities_histogram_left,_1d_densities_histogram_right,
            _1d_densities_histogram
        )
        vpi.histograms.accumulate_radial_densities(
            x[slice_number],
            radial_densities_histogram_maximum_radius,
            radial_densities_histogram
        )
#@-node:gcross.20090819152718.1294:<< Main iteration >>
#@nl

#@<< Compute summary statistics >>
#@+node:gcross.20090819152718.1295:<< Compute summary statistics >>
total_energy_estimates = zeros((2,),dtype=double)
comm.Reduce((energy_estimates,MPI.DOUBLE),(total_energy_estimates,MPI.DOUBLE))

total_left_slice_1d_densities_histogram = compute_total_histogram(left_slice_1d_densities_histogram)
total_center_slice_1d_densities_histogram = compute_total_histogram(center_slice_1d_densities_histogram)
total_right_slice_1d_densities_histogram = compute_total_histogram(right_slice_1d_densities_histogram)

total_left_slice_radial_densities_histogram = compute_total_histogram(left_slice_radial_densities_histogram)
total_center_slice_radial_densities_histogram = compute_total_histogram(center_slice_radial_densities_histogram)
total_right_slice_radial_densities_histogram = compute_total_histogram(right_slice_radial_densities_histogram)

if my_rank == 0:
    #@    << Write results >>
    #@+node:gcross.20090819152718.1296:<< Write results >>
    for move_type_number, move_type_name in enumerate(["Bridge","Rigid"]):
        try:
            print "%s move acceptance rate: %.0f%%" % (
                    move_type_name,
                    float(move_type_accepted_counts[move_type_number])/float(move_type_attempted_counts[move_type_number])*100
                 )
        except ZeroDivisionError:
            print "No %s moves attempted" % (move_type_name.lower())

    print "Observations made:", total_number_of_observations
    print "Energy estimates:", total_energy_estimates/total_number_of_observations

    for slice_name in ["left","center","right"]:
        total_1d_densities_histogram = vars()["total_{slice_name}_slice_1d_densities_histogram".format(**vars())]
        for dimension_number, dimension_name in zip(range(n_dimensions),["x","y","z"]):
            write_histogram(
                "{slice_name}_{dimension_name}_density.dat".format(**vars()),
                total_1d_densities_histogram[dimension_number],
                _1d_densities_histogram_left[dimension_number],
                _1d_densities_histogram_right[dimension_number]
            )
        write_radial_histogram(
            "{slice_name}_radial_density.dat".format(**vars()),
            vars()["total_{slice_name}_slice_radial_densities_histogram".format(**vars())],
            radial_densities_histogram_maximum_radius
        )
    #@-node:gcross.20090819152718.1296:<< Write results >>
    #@nl
#@-node:gcross.20090819152718.1295:<< Compute summary statistics >>
#@nl
#@-node:gcross.20090819152718.1279:@thin run.py
#@-leo
