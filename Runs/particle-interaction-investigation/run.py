#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20090828201103.1780:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20090828201103.1781:<< Imports >>
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
#@-node:gcross.20090828201103.1781:<< Imports >>
#@nl

#@<< System Configuration >>
#@+node:gcross.20090828201103.1782:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_slices": 102,
    "lambda_": 0.5,
    "number_of_dimensions": 2,
    "initial_particle_distribution_size": 4,
    # Angular momentum parameters
    "rotation_plane_axis_1": 1,
    "rotation_plane_axis_2": 2,
    # Run parameters
    "total_number_of_observations": 10000,
    "number_of_prethermalization_steps": 1000,
    # Move parameters
    "dM": 22,
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.05,0.6,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
    # Potential parameters
    "harmonic_oscillator_coefficients": array((1,1),dtype=double),
}
#@-node:gcross.20090828201103.1782:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20090828201103.1783:<< Histogram Configuration >>
_1d_densities_histogram_left  = [-2,-2,-2]
_1d_densities_histogram_right = [+2,+2,+2]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 50
#@-node:gcross.20090828201103.1783:<< Histogram Configuration >>
#@nl
#@<< System properties message >>
#@+node:gcross.20090828201103.2130:<< System properties message >>
system_properties_message = """\
Examining system with
    *) {number_of_particles} particles;
    *) hard core bosons of radius {hard_sphere_radius};
    *) a frame rotating with velocity {frame_angular_velocity};
    *) and {number_of_rotating_particles} particles rotating"""
#@-node:gcross.20090828201103.2130:<< System properties message >>
#@nl

output_root_directory = sys.argv[1]

for hard_sphere_radius in [0,0.1]: 
    for number_of_particles in [10]:
        for frame_angular_velocity in [1]: #8,16]:
            for number_of_rotating_particles in xrange(number_of_particles+1):
                #@                << Run simulation for given parameters >>
                #@+node:gcross.20090828201103.1784:<< Run simulation for given parameters >>
                my_directory = "{output_root_directory}/{number_of_particles}/{hard_sphere_radius}/{frame_angular_velocity}".format(**vars()) 
                if (my_rank == 0):
                    print
                    print system_properties_message.format(**vars())

                for parameter_name in ["hard_sphere_radius","number_of_particles","frame_angular_velocity","number_of_rotating_particles"]:
                    configuration[parameter_name] = vars()[parameter_name]

                system = System(**configuration)
                #@<< Initialize observables >>
                #@+node:gcross.20090828201103.1785:<< Initialize observables >>
                center_slice = system.number_of_slices // 2
                for observable in  [
                    #ParticleSeparationHistogram(center_slice,1,100,"{my_directory}/{number_of_rotating_particles}/particle-separation".format(**vars())),
                    #AverageParticleSeparationEstimate(center_slice,"{my_directory}/particle-separation".format(**vars()),center_slice),
                    TotalEnergyEstimate("{my_directory}/total-energy".format(**vars()),number_of_rotating_particles),
                    ]: system.add_observable(observable)
                #@-node:gcross.20090828201103.1785:<< Initialize observables >>
                #@nl
                system.run()
                system.total_and_write_observables()
                del system.observables
                del system
                gc.collect()
                #@-node:gcross.20090828201103.1784:<< Run simulation for given parameters >>
                #@nl
#@-node:gcross.20090828201103.1780:@thin run.py
#@-leo
