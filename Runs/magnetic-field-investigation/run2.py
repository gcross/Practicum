#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20090922091313.1799:@thin run2.py
#@@first
#@<< Imports >>
#@+node:gcross.20090919132620.2795:<< Imports >>
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
#@-node:gcross.20090919132620.2795:<< Imports >>
#@nl

#@+others
#@+node:gcross.20090919132620.2796:class MagneticFieldPotential
class MagneticFieldPotential(Physics):
    #@    @+others
    #@+node:gcross.20090919132620.2797:hooks
    hooks = ["effective_potentials"] #,"greens_functions"]
    #@-node:gcross.20090919132620.2797:hooks
    #@+node:gcross.20090919132620.2798:__init__
    def __init__(self,system):
        Physics.__init__(self,system)
        self.rotation_plane_axis_1 = system.rotation_plane_axis_1
        self.rotation_plane_axis_2 = system.rotation_plane_axis_2

        self.magnetic_field_strength = system.magnetic_field_strength
        self.rotation_rate = system.rotation_rate

        self.lambda_ = system.lambda_

    #@-node:gcross.20090919132620.2798:__init__
    #@+node:gcross.20090919132620.2799:accumulate_potential
    def accumulate_potential(self,x,xij2,U,gradU):

        numeric_gradients = zeros(x.shape,dtype=double,order='Fortran')

        vpif.angular_momentum.accumulate_gradient_feynman(
            x,
            self.rotation_rate,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            numeric_gradients
        )

        vpif.angular_momentum.accumulate_magnetic_field_phase(
            x,
            self.magnetic_field_strength,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            numeric_gradients
        )

        vpif.angular_momentum.accumulate_effective_potential(
            numeric_gradients,
            self.lambda_,
            U
        )
    #@-node:gcross.20090919132620.2799:accumulate_potential
    #@+node:gcross.20090919132620.2800:compute_greens_function
    def compute_greens_function(self,
            x,xij2,
            U,gradU2,
            lam,dt,
            slice_start,slice_end,
            particle_number
        ): return \
            vpif.angular_momentum.compute_greens_function(
                x,
                lam, dt,
                slice_start, slice_end,
                self.number_of_rotating_particles,
                self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            )
    #@-node:gcross.20090919132620.2800:compute_greens_function
    #@-others
#@-node:gcross.20090919132620.2796:class MagneticFieldPotential
#@-others

#@<< System Configuration >>
#@+node:gcross.20090919132620.2801:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_slices": 202,
    "lambda_": 0.5,
    "number_of_dimensions": 2,
    "initial_particle_distribution_size": 1,
    # Angular momentum parameters
    "rotation_plane_axis_1": 1,
    "rotation_plane_axis_2": 2,
    # Run parameters
    "total_number_of_observations": 100000,
    "number_of_prethermalization_steps": 2000,
    # Move parameters
    "dM": 96,
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.01,0.025,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
}
configuration["harmonic_oscillator_coefficients"] = \
    array([1.0,]*configuration["number_of_dimensions"],dtype=double)
#@nonl
#@-node:gcross.20090919132620.2801:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20090919132620.2802:<< Histogram Configuration >>
_1d_densities_histogram_left  = [-2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_right = [+2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 50
#@-node:gcross.20090919132620.2802:<< Histogram Configuration >>
#@nl
#@<< System properties message >>
#@+node:gcross.20090919132620.2803:<< System properties message >>
system_properties_message = """\
Examining system with
    *) {number_of_particles} particles;
    *) a magnetic field with strength {magnetic_field_strength};
    *) and a feynman phase multipier of {rotation_rate}"""
#@-node:gcross.20090919132620.2803:<< System properties message >>
#@nl

output_root_directory = sys.argv[1]

for number_of_particles in [10]:
    for magnetic_field_strength in [-4,-8]:
        for rotation_rate in xrange(0,201+1,40):
            #@            << Run simulation for given parameters >>
            #@+node:gcross.20090919132620.2804:<< Run simulation for given parameters >>
            my_directory = "{output_root_directory}/{number_of_particles}/{magnetic_field_strength}".format(**vars()) 
            if (my_rank == 0):
                print
                print system_properties_message.format(**vars())

            for parameter_name in ["number_of_particles","magnetic_field_strength","rotation_rate"]:
                configuration[parameter_name] = vars()[parameter_name]

            system = System(**configuration)
            #@<< Initialize physics >>
            #@+node:gcross.20090919132620.2805:<< Initialize physics >>
            for physics in [
                HarmonicOscillator,
                MagneticFieldPotential,
                SecondOrderGreensFunction,
                ]: system.add_physics(physics)
            #@-node:gcross.20090919132620.2805:<< Initialize physics >>
            #@nl
            #@<< Initialize observables >>
            #@+node:gcross.20090919132620.2806:<< Initialize observables >>
            for slice_name, slice_number in [("left",0),("center",system.center_slice_number),("right",system.number_of_slices-1)]:
                density_slice_subdirectory = "{my_directory}/{rotation_rate}/{slice_name}".format(**vars())
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
                            rotation_rate
                        ),
                        PhysicalPotentialSliceEnergyEstimate(
                            slice_number,
                            "{my_directory}/{slice_name}/physical-potential".format(**vars()),
                            rotation_rate
                        ),
                        TotalPotentialSliceEnergyEstimate(
                            slice_number,
                            "{my_directory}/{slice_name}/total-potential".format(**vars()),
                            rotation_rate
                        ),
                    ] + [
                        AverageAxialDistanceEstimate(
                            axis_number,
                            slice_number,
                            "{my_directory}/{slice_name}/average-{axis_name}".format(**vars()),
                            rotation_rate
                        ) for (axis_number,axis_name) in enumerate(["x","y","x","w"][:system.number_of_dimensions])
                    ] + [
                        AverageRadiusEstimate(
                            slice_number,
                            "{my_directory}/{slice_name}/average-radius".format(**vars()),
                            rotation_rate
                        ),
                    ]:
                    system.add_observable(observable)

            center_slice = system.number_of_slices // 2
            for observable in  [
                TotalEnergyEstimate("{my_directory}/total-energy".format(**vars()),rotation_rate)
                ]: system.add_observable(observable)
            #@-node:gcross.20090919132620.2806:<< Initialize observables >>
            #@nl
            system.run()
            system.total_and_write_observables()
            del system.observables
            del system
            gc.collect()
            #@-node:gcross.20090919132620.2804:<< Run simulation for given parameters >>
            #@nl
#@-node:gcross.20090922091313.1799:@thin run2.py
#@-leo
