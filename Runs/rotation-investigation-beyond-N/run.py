#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20090918204102.1765:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20090918204102.1766:<< Imports >>
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

from numpy import log
#@-node:gcross.20090918204102.1766:<< Imports >>
#@nl

#@+others
#@+node:gcross.20090918204102.1767:class RotationEffectivePotential
class RotationEffectivePotential(Physics):
    #@    @+others
    #@+node:gcross.20090918204102.1768:hooks
    hooks = ["effective_potentials"] #,"greens_functions"]
    #@-node:gcross.20090918204102.1768:hooks
    #@+node:gcross.20090918204102.1769:__init__
    def __init__(self,system):
        Physics.__init__(self,system)
        self.rotation_plane_axis_1 = system.rotation_plane_axis_1
        self.rotation_plane_axis_2 = system.rotation_plane_axis_2
        self.frame_angular_velocity = system.frame_angular_velocity

        self.number_of_rotating_particles = system.number_of_rotating_particles

        self.lambda_ = system.lambda_

    #@-node:gcross.20090918204102.1769:__init__
    #@+node:gcross.20090918204102.1770:accumulate_potential
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
    #@-node:gcross.20090918204102.1770:accumulate_potential
    #@+node:gcross.20090918204102.1771:compute_greens_function
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
    #@-node:gcross.20090918204102.1771:compute_greens_function
    #@-others
#@-node:gcross.20090918204102.1767:class RotationEffectivePotential
#@+node:gcross.20090918204102.1772:class LogDistanceToWallHistogram
class LogDistanceToWallHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090918204102.1773:(fields)
    left = 0
    #@-node:gcross.20090918204102.1773:(fields)
    #@+node:gcross.20090918204102.1774:__init__
    def __init__(self,slice_number,maximum_distance,number_of_bins,filename):
        self.right = maximum_distance
        self.number_of_bins = number_of_bins
        self.dndx = float(number_of_bins) / maximum_distance
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090918204102.1774:__init__
    #@+node:gcross.20090918204102.1775:update
    def update(self):
        system = self.system
        distance = \
            vpif.angular_momentum.estimate_distance_to_node(
                system.x[self.slice_number],
                system.number_of_rotating_particles,
                system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            )/sqrt(system.lambda_*system.move_type_differentials[0])
        vpif.histograms.place_in_bin(
            distance,
            0, # offset
            self.dndx,
            self.histogram
        )
    #@-node:gcross.20090918204102.1775:update
    #@-others
#@-node:gcross.20090918204102.1772:class LogDistanceToWallHistogram
#@+node:gcross.20090918204102.1776:class GreenCorrectionHistogram
class GreenCorrectionHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090918204102.1777:(fields)
    left = 0
    #@-node:gcross.20090918204102.1777:(fields)
    #@+node:gcross.20090918204102.1778:__init__
    def __init__(self,slice_number,maximum_distance,number_of_bins,filename):
        self.right = maximum_distance
        self.number_of_bins = number_of_bins
        self.dndx = float(number_of_bins) / maximum_distance
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090918204102.1778:__init__
    #@+node:gcross.20090918204102.1779:update
    def update(self):
        system = self.system
        value = \
            vpif.angular_momentum.compute_greens_function(
                system.x,
                system.lambda_, system.move_type_differentials[0],
                1, system.number_of_slices,
                system.number_of_rotating_particles,
                system.rotation_plane_axis_1,system.rotation_plane_axis_2,
            )
        vpif.histograms.place_in_bin(
            exp(value),
            0, # offset
            self.dndx,
            self.histogram
        )
    #@-node:gcross.20090918204102.1779:update
    #@-others
#@-node:gcross.20090918204102.1776:class GreenCorrectionHistogram
#@-others

#@<< System Configuration >>
#@+node:gcross.20090918204102.1780:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_slices": 202,
    "lambda_": 0.5,
    "number_of_dimensions": 3,
    "initial_particle_distribution_size": 1,
    # Angular momentum parameters
    "rotation_plane_axis_1": 1,
    "rotation_plane_axis_2": 2,
    # Run parameters
    "total_number_of_observations": 100000,
    "number_of_prethermalization_steps": 4000,
    # Move parameters
    "dM": 96,
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.01,0.3,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
}
configuration["harmonic_oscillator_coefficients"] = \
    array([1.0,]*configuration["number_of_dimensions"],dtype=double)
#@nonl
#@-node:gcross.20090918204102.1780:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20090918204102.1781:<< Histogram Configuration >>
_1d_densities_histogram_left  = [-2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_right = [+2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 50
#@-node:gcross.20090918204102.1781:<< Histogram Configuration >>
#@nl
#@<< System properties message >>
#@+node:gcross.20090918204102.1782:<< System properties message >>
system_properties_message = """\
Examining system with
    *) {number_of_particles} particles;
    *) a frame rotating with velocity {frame_angular_velocity};
    *) and {number_of_rotating_particles} particles rotating"""
#@-node:gcross.20090918204102.1782:<< System properties message >>
#@nl

output_root_directory = sys.argv[1]

for number_of_particles in [6]:
    for frame_angular_velocity in [1,2]:
        for number_of_rotating_particles in xrange(0,20*number_of_particles+1):
            #@            << Run simulation for given parameters >>
            #@+node:gcross.20090918204102.1783:<< Run simulation for given parameters >>
            my_directory = "{output_root_directory}/{number_of_particles}/{frame_angular_velocity}".format(**vars()) 
            if (my_rank == 0):
                print
                print system_properties_message.format(**vars())

            for parameter_name in ["number_of_particles","frame_angular_velocity","number_of_rotating_particles"]:
                configuration[parameter_name] = vars()[parameter_name]

            system = System(**configuration)
            #@<< Initialize physics >>
            #@+node:gcross.20090918204102.1784:<< Initialize physics >>
            for physics in [
                HarmonicOscillator,
                RotationEffectivePotential,
                SecondOrderGreensFunction,
                ]: system.add_physics(physics)
            #@-node:gcross.20090918204102.1784:<< Initialize physics >>
            #@nl
            #@<< Initialize observables >>
            #@+node:gcross.20090918204102.1785:<< Initialize observables >>
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
            #@-node:gcross.20090918204102.1785:<< Initialize observables >>
            #@nl
            system.run()
            system.total_and_write_observables()
            del system.observables
            del system
            gc.collect()
            #@-node:gcross.20090918204102.1783:<< Run simulation for given parameters >>
            #@nl
#@-node:gcross.20090918204102.1765:@thin run.py
#@-leo
