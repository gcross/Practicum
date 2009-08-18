#@+leo-ver=4-thin
#@+node:gcross.20090818114910.1259:@thin probe.py
import sys
sys.path.append("lib")

#@<< Imports >>
#@+node:gcross.20090818114910.1260:<< Imports >>
import vpi

from numpy import *
from numpy.random import rand

from mpi4py import MPI

comm = MPI.COMM_WORLD

number_of_processors = comm.Get_size()
#@-node:gcross.20090818114910.1260:<< Imports >>
#@nl

#@<< Configuration >>
#@+node:gcross.20090818114910.1261:<< Configuration >>
# System parameters
n_particles = 4
n_slices = 54
lambda_ = 0.5

# Run parameters
total_number_of_observations = 10000
number_of_thermalization_steps = 1000

# Angular momentum parameters
rotation_plane_axis_1 = 1
rotation_plane_axis_2 = 2

# Move parameters
dM = 22
move_type_probabilities = [0.9,0.1,0]
move_type_differentials = [0.1,0.6,0]
low_swap_dimension = 1
high_swap_dimension= 3

# Potential parameters
harmonic_oscillator_coefficients = (1,1,1)

# Miscellaneous
#@-node:gcross.20090818114910.1261:<< Configuration >>
#@nl

#@+others
#@+node:gcross.20090818114910.1274:probe
def probe(filename,frame_angular_velocity,N_rotating_particles):
    #@    << Initialization >>
    #@+node:gcross.20090818114910.1268:<< Initialization >>
    n_dimensions = 3

    assert (n_slices % 2 == 0 and n_slices % 4 == 2)
    c_slice = n_slices / 2

    x = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')
    x[...] = rand(1,n_particles,n_dimensions)

    xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
    vpi.xij.update_xij(xij2,x)

    U = zeros((n_slices,n_particles),dtype=double,order='Fortran')
    gradU2 = zeros((n_slices),dtype=double,order='Fortran')

    slice_move_attempted_counts = zeros((n_slices,),'i')
    slice_move_accepted_counts = zeros((n_slices,),'i')

    move_type_attempted_counts = zeros((3,),'i')
    move_type_accepted_counts = zeros((3,),'i')

    U_weights, gU2_weights = vpi.gfn.initialize_4th_order_weights(n_slices)

    n_trials = n_particles*n_slices/dM
    use_4th_order_green_function = False
    #@-node:gcross.20090818114910.1268:<< Initialization >>
    #@nl

    #@    @+others
    #@+node:gcross.20090818114910.1263:Functions
    #@+node:gcross.20090818114910.1264:Potential function
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
    #@-node:gcross.20090818114910.1264:Potential function
    #@+node:gcross.20090818114910.1265:Trial function
    def trial_function(x,xij2):
        return -sum(dot(x**2,harmonic_oscillator_coefficients))/2
    #@-node:gcross.20090818114910.1265:Trial function
    #@+node:gcross.20090818114910.1266:Trial derivatives
    def trial_derivatives(x,xij2):
        gradient_of_log_trial_fn = -x*harmonic_oscillator_coefficients
        laplacian_of_log_trial_fn = -sum(harmonic_oscillator_coefficients)*x.shape[0]
        return gradient_of_log_trial_fn, laplacian_of_log_trial_fn
    #@-node:gcross.20090818114910.1266:Trial derivatives
    #@+node:gcross.20090818114910.1267:Observables
    def compute_energy_estimates():
        estimates = []
        for i in [0,-1]:
            gradient_of_log_trial_fn, laplacian_of_log_trial_fn = trial_derivatives(x[i],xij2[i])
            estimates.append(
                vpi.observables.compute_local_energy_estimate(
                    U[0],
                    gradient_of_log_trial_fn, laplacian_of_log_trial_fn,
                    lambda_,
                )
            )
        return estimates
    #@-node:gcross.20090818114910.1267:Observables
    #@-node:gcross.20090818114910.1263:Functions
    #@-others

    #@    << Main iteration >>
    #@+node:gcross.20090818114910.1269:<< Main Iteration >>
    vpi.thermalize.thermalize_path(
        x,xij2,
        U,gradU2,
        n_trials*number_of_thermalization_steps,
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

    number_of_observations = total_number_of_observations // number_of_processors

    for _ in xrange(number_of_observations):
        vpi.thermalize.thermalize_path(
            x,xij2,
            U,gradU2,
            n_trials,
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

    #@-node:gcross.20090818114910.1269:<< Main Iteration >>
    #@nl

    #@    << Write out result >>
    #@+node:gcross.20090818114910.1271:<< Write out result >>
    total_energy_estimates = zeros((2,),dtype=double)

    comm.Reduce((energy_estimates,MPI.DOUBLE),(total_energy_estimates,MPI.DOUBLE))

    if comm.Get_rank() == 0:
        total_energy_estimates /= (comm.Get_size()*number_of_observations)
        total_energy_estimate = sum(total_energy_estimates)/len(total_energy_estimates)
        total_energy_estimate_difference = abs(total_energy_estimates[0]-total_energy_estimates[1])
        print "%.0f%%; %f,%i: E=%f, dE=%f" % (
                float(move_type_accepted_counts[0])/float(move_type_attempted_counts[0])*100,
                frame_angular_velocity,
                N_rotating_particles,
                total_energy_estimate,
                total_energy_estimate_difference,
            )
        with open("{0}-{1}-{2}.dat".format(filename,frame_angular_velocity,N_rotating_particles),"w") as f:
            print >> f, " ".join(map(str,[frame_angular_velocity,N_rotating_particles,total_energy_estimate]))
    #@-node:gcross.20090818114910.1271:<< Write out result >>
    #@nl
#@-node:gcross.20090818114910.1274:probe
#@-others

if __name__ == "__main__":
    filename = sys.argv[1]
    frame_angular_velocity = float(sys.argv[2])
    N_rotating_particles = int(sys.argv[3])
    probe(filename,frame_angular_velocity,N_rotating_particles)
#@-node:gcross.20090818114910.1259:@thin probe.py
#@-leo
