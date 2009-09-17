#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20090911091023.2427:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20090911091023.2428:<< Imports >>
import gc

import sys
sys.path.append("lib")

try:
    import psyco
    psyco.full()
except ImportError:
    pass

from vpi import *
import vpif

import itertools
#@-node:gcross.20090911091023.2428:<< Imports >>
#@nl

#@+others
#@+node:gcross.20090911091023.2443:class RotationEffectivePotential
class RotationEffectivePotential(Physics):
    #@    @+others
    #@+node:gcross.20090911091023.2444:hooks
    hooks = ["effective_potentials"]
    #@-node:gcross.20090911091023.2444:hooks
    #@+node:gcross.20090911091023.2445:__init__
    def __init__(self,system):
        Physics.__init__(self,system)
        self.rotation_plane_axis_1 = system.rotation_plane_axis_1
        self.rotation_plane_axis_2 = system.rotation_plane_axis_2
        self.frame_angular_velocity = system.frame_angular_velocity

        self.number_of_rotating_particles = system.number_of_rotating_particles

        self.lambda_ = system.lambda_

    #@-node:gcross.20090911091023.2445:__init__
    #@+node:gcross.20090911091023.2446:accumulate_potential
    def accumulate_potential(self,x,xij2,U,gradU):

        numeric_gradients = zeros(x.shape,dtype=double,order='Fortran')

        vpif.angular_momentum.accumulate_gradient_fancy(
            x,
            self.number_of_rotating_particles,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            numeric_gradients
        )

        vpif.angular_momentum.accumulate_effective_potential(
            x, numeric_gradients,
            self.frame_angular_velocity,self.lambda_,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            U
        )
    #@-node:gcross.20090911091023.2446:accumulate_potential
    #@-others
#@-node:gcross.20090911091023.2443:class RotationEffectivePotential
#@-others

#@<< System Configuration >>
#@+node:gcross.20090911091023.2433:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_slices": 102,
    "lambda_": 0.5,
    "number_of_dimensions": 3,
    "initial_particle_distribution_size": 1,
    # Angular momentum parameters
    "rotation_plane_axis_1": 1,
    "rotation_plane_axis_2": 2,
    # Run parameters
    "total_number_of_observations": 10000,
    "number_of_prethermalization_steps": 2000,
    # Move parameters
    "dM": 16,
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.05,0.6,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
}
configuration["harmonic_oscillator_coefficients"] = \
    array([1.0,]*configuration["number_of_dimensions"],dtype=double)
#@nonl
#@-node:gcross.20090911091023.2433:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20090911091023.2434:<< Histogram Configuration >>
_1d_densities_histogram_left  = [-2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_right = [+2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 50
#@-node:gcross.20090911091023.2434:<< Histogram Configuration >>
#@nl
#@<< System properties message >>
#@+node:gcross.20090911091023.2435:<< System properties message >>
system_properties_message = """\
Examining system with
    *) {number_of_particles} particles;
    *) a frame rotating with velocity {frame_angular_velocity};
    *) and {number_of_rotating_particles} particles rotating"""
#@-node:gcross.20090911091023.2435:<< System properties message >>
#@nl

output_root_directory = sys.argv[1]

for number_of_particles in [6]:
    for frame_angular_velocity in [0.5,1]: #[1]: #8,16]:
        for number_of_rotating_particles in xrange(number_of_particles+1):
            #@            << Run simulation for given parameters >>
            #@+node:gcross.20090911091023.2436:<< Run simulation for given parameters >>
            my_directory = "{output_root_directory}/{number_of_particles}/{frame_angular_velocity}".format(**vars()) 
            if (my_rank == 0):
                print
                print system_properties_message.format(**vars())

            for parameter_name in ["number_of_particles","frame_angular_velocity","number_of_rotating_particles"]:
                configuration[parameter_name] = vars()[parameter_name]

            system = System(**configuration)
            #@<< Initialize physics >>
            #@+node:gcross.20090911091023.2437:<< Initialize physics >>
            for physics in [
                HarmonicOscillator,
                RotationEffectivePotential,
                SecondOrderGreensFunction,
                ]: system.add_physics(physics)
            #@-node:gcross.20090911091023.2437:<< Initialize physics >>
            #@nl
            #@<< Initialize observables >>
            #@+node:gcross.20090911091023.2438:<< Initialize observables >>
            for slice_name, slice_number in [("left",0),("center",system.center_slice_number),("right",system.number_of_slices-1)]:
                density_slice_subdirectory = "{my_directory}/{number_of_rotating_particles}/{slice_name}".format(**vars())
                for observable in [
                        PositionDensity1DHistogram(
                            slice_number,
                            _1d_densities_histogram_left,
                            _1d_densities_histogram_right,
                            _1d_densities_histogram_bin_count,
                            [density_slice_subdirectory + "/1d-densities/" + label for label in ["x","y","z","w"][:system.number_of_dimensions]]
                        ),
                        RadialDensityHistogram(
                            slice_number,
                            radial_densities_histogram_maximum_radius,
                            radial_densities_histogram_bin_count,
                            density_slice_subdirectory + "/radial-density"
                        ),
                        EffectivePotentialSliceEnergyEstimate(
                            slice_number,
                            "{my_directory}/{slice_name}/effective-potential".format(**vars()),
                            number_of_rotating_particles
                        ),
                        PhysicalPotentialSliceEnergyEstimate(
                            slice_number,
                            "{my_directory}/{slice_name}/physical-potential".format(**vars()),
                            number_of_rotating_particles
                        ),
                        TotalPotentialSliceEnergyEstimate(
                            slice_number,
                            "{my_directory}/{slice_name}/total-potential".format(**vars()),
                            number_of_rotating_particles
                        ),
                    ] + [
                        AverageAxialDistanceEstimate(
                            axis_number,
                            slice_number,
                            "{my_directory}/{slice_name}/average-{axis_name}".format(**vars()),
                            number_of_rotating_particles
                        ) for (axis_number,axis_name) in enumerate(["x","y","x","w"][:system.number_of_dimensions])
                    ] + [
                        AverageRadiusEstimate(
                            slice_number,
                            "{my_directory}/{slice_name}/average-radius".format(**vars()),
                            number_of_rotating_particles
                        ),
                    ]:
                    system.add_observable(observable)

            center_slice = system.number_of_slices // 2
            for observable in  [
                #ParticleSeparationHistogram(center_slice,1,100,"{my_directory}/{number_of_rotating_particles}/particle-separation".format(**vars())),
                #AverageParticleSeparationEstimate(center_slice,"{my_directory}/particle-separation".format(**vars()),center_slice),
                TotalEnergyEstimate("{my_directory}/total-energy".format(**vars()),number_of_rotating_particles),
                ]: system.add_observable(observable)
            #@-node:gcross.20090911091023.2438:<< Initialize observables >>
            #@nl
            system.run()
            system.total_and_write_observables()
            del system.observables
            del system
            gc.collect()
            #@-node:gcross.20090911091023.2436:<< Run simulation for given parameters >>
            #@nl
#@-node:gcross.20090911091023.2427:@thin run.py
#@-leo
