#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20090827130017.1769:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20090827130017.1770:<< Imports >>
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
#@-node:gcross.20090827130017.1770:<< Imports >>
#@nl

#@<< System Configuration >>
#@+node:gcross.20090827130017.1771:<< System Configuration >>
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
#@-node:gcross.20090827130017.1771:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20090827130017.1772:<< Histogram Configuration >>
_1d_densities_histogram_left  = [-2,-2,-2]
_1d_densities_histogram_right = [+2,+2,+2]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 50
#@-node:gcross.20090827130017.1772:<< Histogram Configuration >>
#@nl

output_root_directory = sys.argv[1]

for number_of_particles in [6]:
    for frame_angular_velocity in [4,8,16]: #2,4,8,16]:
        for number_of_rotating_particles in xrange(number_of_particles+1):
            #@            << Run simulation for given parameters >>
            #@+node:gcross.20090827130017.1773:<< Run simulation for given parameters >>
            my_directory = "{output_root_directory}/{number_of_particles}/{frame_angular_velocity}".format(**vars()) 
            if (my_rank == 0):
                print
                print "Examining system with {number_of_particles} particles in a frame rotating with velocity {frame_angular_velocity} and with {number_of_rotating_particles} particles rotating".format(**vars())

            for parameter_name in ["number_of_particles","frame_angular_velocity","number_of_rotating_particles"]:
                configuration[parameter_name] = vars()[parameter_name]

            system = System(**configuration)
            #@<< Initialize observables >>
            #@+node:gcross.20090827130017.1774:<< Initialize observables >>
            for observable in  [
                TotalEnergyEstimate("{my_directory}/total-energy".format(**vars()),number_of_rotating_particles),
                ]: system.add_observable(observable)
            #@-node:gcross.20090827130017.1774:<< Initialize observables >>
            #@nl
            system.run()
            system.total_and_write_observables()
            del system.observables
            del system
            gc.collect()
            #@-node:gcross.20090827130017.1773:<< Run simulation for given parameters >>
            #@nl
#@-node:gcross.20090827130017.1769:@thin run.py
#@-leo
