#@+leo-ver=4-thin
#@+node:gcross.20090819152718.1315:@thin probe.py
import sys
sys.path.append("../../lib")

#@<< Imports >>
#@+node:gcross.20090819152718.1316:<< Imports >>
import vpi

from numpy import *
from numpy.random import rand

#@-node:gcross.20090819152718.1316:<< Imports >>
#@nl

#@<< Configuration >>
#@+node:gcross.20090819152718.1317:<< Configuration >>
# System parameters
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
use_4th_order_green_function = False
#@-node:gcross.20090819152718.1317:<< Configuration >>
#@nl

#@<< MPI Initialization >>
#@+node:gcross.20090819152718.1334:<< MPI Initialization >>
from mpi4py import MPI

comm = MPI.COMM_WORLD

number_of_processors = comm.Get_size()
my_rank = comm.Get_rank()

number_of_observations = total_number_of_observations // number_of_processors + 1

total_number_of_observations = number_of_observations * number_of_processors
#@nonl
#@-node:gcross.20090819152718.1334:<< MPI Initialization >>
#@nl

#@+others
#@+node:gcross.20090819152718.1333:probe
def probe(n_particles, frame_angular_velocity,N_rotating_particles):
    #@    @+others
    #@+node:gcross.20090819152718.1318:Functions
    #@+node:gcross.20090819152718.1319:Physics
    #@+node:gcross.20090819152718.1320:Potential function
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
    #@-node:gcross.20090819152718.1320:Potential function
    #@+node:gcross.20090819152718.1321:Trial function
    def trial_function(x,_):
        return -sum(dot(x**2,harmonic_oscillator_coefficients))/2
    #@-node:gcross.20090819152718.1321:Trial function
    #@+node:gcross.20090819152718.1322:Trial derivatives
    def trial_derivatives(x,xij2):
        gradient_of_log_trial_fn = -x*harmonic_oscillator_coefficients
        laplacian_of_log_trial_fn = -sum(harmonic_oscillator_coefficients)*x.shape[0]
        return gradient_of_log_trial_fn, laplacian_of_log_trial_fn
    #@-node:gcross.20090819152718.1322:Trial derivatives
    #@-node:gcross.20090819152718.1319:Physics
    #@+node:gcross.20090819152718.1323:Observables
    #@+node:gcross.20090819152718.1324:Energy
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
    #@-node:gcross.20090819152718.1324:Energy
    #@-node:gcross.20090819152718.1323:Observables
    #@-node:gcross.20090819152718.1318:Functions
    #@-others

    #@    << Initialization >>
    #@+node:gcross.20090819152718.1329:<< Initialization >>
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
    #@-node:gcross.20090819152718.1329:<< Initialization >>
    #@nl

    #@    << Main iteration >>
    #@+node:gcross.20090819152718.1330:<< Main iteration >>
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
    #@-node:gcross.20090819152718.1330:<< Main iteration >>
    #@nl

    #@    << Compute summary statistics >>
    #@+node:gcross.20090819152718.1331:<< Compute summary statistics >>
    total_energy_estimates = zeros((2,),dtype=double)
    comm.Reduce((energy_estimates,MPI.DOUBLE),(total_energy_estimates,MPI.DOUBLE))

    total_energy_estimates /= total_number_of_observations

    if(my_rank == 0):
        return mean(total_energy_estimates), std(total_energy_estimates)
    #@-node:gcross.20090819152718.1331:<< Compute summary statistics >>
    #@nl
#@-node:gcross.20090819152718.1333:probe
#@-others
#@-node:gcross.20090819152718.1315:@thin probe.py
#@-leo
