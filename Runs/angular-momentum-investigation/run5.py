#@+leo-ver=4-thin
#@+node:gcross.20090826111206.1451:@thin run5.py
#@<< Imports >>
#@+node:gcross.20090826111206.1452:<< Imports >>
import gc

import sys
sys.path.append("lib")

try:
    import psyco
    psyco.full()
except ImportError:
    pass

from system import *

import itertools
#@-node:gcross.20090826111206.1452:<< Imports >>
#@nl

#@<< System Configuration >>
#@+node:gcross.20090826111206.1453:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_slices": 102,
    "lambda_": 0.5,
    "number_of_dimensions": 3,
    # Angular momentum parameters
    "rotation_plane_axis_1": 1,
    "rotation_plane_axis_2": 2,
    # Run parameters
    "total_number_of_observations": 100000,
    "number_of_prethermalization_steps": 1000,
    # Move parameters
    "dM": 22,
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.1,0.6,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
    # Potential parameters
    "harmonic_oscillator_coefficients": array((1,1,1),dtype=double),
    # Miscellaneous
    "use_4th_order_green_function": False,
}
#@-node:gcross.20090826111206.1453:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20090826111206.1454:<< Histogram Configuration >>
_1d_densities_histogram_left  = [-2,-2,-2]
_1d_densities_histogram_right = [+2,+2,+2]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 50
#@-node:gcross.20090826111206.1454:<< Histogram Configuration >>
#@nl

output_root_directory = sys.argv[1]

for number_of_particles in [6]:
    for frame_angular_velocity in [0.4]:
        for number_of_rotating_particles in xrange(number_of_particles+1):
            #@            << Run simulation for given parameters >>
            #@+node:gcross.20090826111206.1455:<< Run simulation for given parameters >>
            my_directory = "{output_root_directory}/{number_of_particles}/{frame_angular_velocity}".format(**vars()) 
            if (my_rank == 0):
                print
                print "Examining system with {number_of_particles} particles in a frame rotating with velocity {frame_angular_velocity} and with {number_of_rotating_particles} particles rotating".format(**vars())

            for parameter_name in ["number_of_particles","frame_angular_velocity","number_of_rotating_particles"]:
                configuration[parameter_name] = vars()[parameter_name]

            system = System(**configuration)
            #@<< Initialize observables >>
            #@+node:gcross.20090826111206.1456:<< Initialize observables >>
            for slice_name, slice_number in [("center",system.center_slice_number)]:
                density_slice_subdirectory = "{my_directory}/{number_of_rotating_particles}/{slice_name}".format(**vars())
                for observable in [
                        AngularVelocitySquaredHistogram(
                            slice_number,
                            0,1,
                            50,
                            density_slice_subdirectory + "/angular-velocities-squared"
                        ),
                        AverageAngularVelocitySquaredEstimate(
                            slice_number,
                            "{my_directory}/{slice_name}/average-angular-velocity-squared".format(**vars()),
                            number_of_rotating_particles
                        ),
                    ]:
                    system.add_observable(observable)

            for observable in  [
                TotalEnergyEstimate("{my_directory}/total-energy".format(**vars()),number_of_rotating_particles),
                ]: system.add_observable(observable)
            #@-node:gcross.20090826111206.1456:<< Initialize observables >>
            #@nl
            system.run()
            system.total_and_write_observables()
            del system.observables
            del system
            gc.collect()
            #@-node:gcross.20090826111206.1455:<< Run simulation for given parameters >>
            #@nl
#@-node:gcross.20090826111206.1451:@thin run5.py
#@-leo
