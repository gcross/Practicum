#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20090908092630.1851:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20090908092630.1852:<< Imports >>
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
#@-node:gcross.20090908092630.1852:<< Imports >>
#@nl

#@+others
#@+node:gcross.20090908092630.1853:class FeynmanEffectivePotentialEstimate
class FeynmanEffectivePotentialEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090908092630.1854:update
    def update(self):
        system = self.system

        x = system.x[self.slice_number:self.slice_number+1]

        numeric_gradients = zeros(x.shape,dtype=double,order='Fortran')
        vpif.angular_momentum.accumulate_gradient_feynman(
            x,
            system.rotation_rate,
            system.rotation_plane_axis_1,system.rotation_plane_axis_2,
            numeric_gradients
        )
        U = zeros(x.shape[:2],dtype=double,order='Fortran')
        vpif.angular_momentum.accumulate_effective_potential(
            numeric_gradients,
            system.frame_angular_velocity,system.lambda_,
            U
        )

        self.add(sum(U))
    #@-node:gcross.20090908092630.1854:update
    #@-others
#@-node:gcross.20090908092630.1853:class FeynmanEffectivePotentialEstimate
#@+node:gcross.20090908092630.1855:class BackflowEffectivePotentialEstimate
class BackflowEffectivePotentialEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090908092630.1856:update
    def update(self):
        system = self.system

        x = system.x[self.slice_number:self.slice_number+1]
        xij2 = system.xij2[self.slice_number:self.slice_number+1]

        backflow_coefficient = -system.backflow_coefficient*system.rotation_rate
        numeric_gradients = vpif.hard_sphere_interaction.compute_gradient_backflow(
            x,xij2,
            system.hard_sphere_radius,backflow_coefficient,
            system.rotation_plane_axis_1,system.rotation_plane_axis_2
        )

        U = zeros(x.shape[:2],dtype=double,order='Fortran')
        vpif.angular_momentum.accumulate_effective_potential(
            numeric_gradients,
            system.frame_angular_velocity,system.lambda_,
            U
        )

        self.add(sum(U))
    #@-node:gcross.20090908092630.1856:update
    #@-others
#@-node:gcross.20090908092630.1855:class BackflowEffectivePotentialEstimate
#@+node:gcross.20090908092630.1857:class RotationEffectivePotential
class RotationEffectivePotential(Physics):
    #@    @+others
    #@+node:gcross.20090908092630.1858:hooks
    hooks = ["effective_potentials"]
    #@-node:gcross.20090908092630.1858:hooks
    #@+node:gcross.20090908092630.1859:__init__
    def __init__(self,system):
        Physics.__init__(self,system)
        self.rotation_plane_axis_1 = system.rotation_plane_axis_1
        self.rotation_plane_axis_2 = system.rotation_plane_axis_2
        self.frame_angular_velocity = system.frame_angular_velocity

        self.rotation_rate = float(system.number_of_rotating_particles) / system.number_of_particles
        self.backflow_coefficient = -system.backflow_coefficient*self.rotation_rate

        self.lambda_ = system.lambda_

        self.hard_sphere_radius = system.hard_sphere_radius
    #@-node:gcross.20090908092630.1859:__init__
    #@+node:gcross.20090908092630.1860:accumulate_potential
    def accumulate_potential(self,x,xij2,U,gradU):

        numeric_gradients = vpif.hard_sphere_interaction.compute_gradient_backflow(
            x,xij2,
            self.hard_sphere_radius,self.backflow_coefficient,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2
        )

    #    numeric_gradients = zeros(x.shape,dtype=double,order='Fortran')

        vpif.angular_momentum.accumulate_gradient_feynman(
            x,
            self.rotation_rate,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            numeric_gradients
        )
        vpif.angular_momentum.accumulate_effective_potential(
            numeric_gradients,
            self.frame_angular_velocity,self.lambda_,
            U
        )
    #@-node:gcross.20090908092630.1860:accumulate_potential
    #@-others
#@-node:gcross.20090908092630.1857:class RotationEffectivePotential
#@-others

#@<< System Configuration >>
#@+node:gcross.20090908092630.1861:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_slices": 142,
    "lambda_": 0.5,
    "number_of_dimensions": 3,
    "initial_particle_distribution_size": 8,
    # Angular momentum parameters
    "rotation_plane_axis_1": 1,
    "rotation_plane_axis_2": 2,
    "frame_angular_velocity": 0,
    # Run parameters
    "total_number_of_observations": 2, #40000,
    "number_of_prethermalization_steps": 10, #2000,
    # Move parameters
    "dM": 24,
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.03,0.6,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
}


#@-node:gcross.20090908092630.1861:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20090908092630.1862:<< Histogram Configuration >>
_1d_densities_histogram_left  = [-2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_right = [+2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 50
#@-node:gcross.20090908092630.1862:<< Histogram Configuration >>
#@nl
#@<< System properties message >>
#@+node:gcross.20090908092630.1863:<< System properties message >>
system_properties_message = """\
Examining system with
    *) {number_of_particles} particles;
    *) harmonic oscillator strength {harmonic_oscillator_strength};
    *) hard core bosons of radius {hard_sphere_radius};
    *) and a backflow correction coefficient of {backflow_coefficient};
"""
#@nonl
#@-node:gcross.20090908092630.1863:<< System properties message >>
#@nl

output_root_directory = sys.argv[1]

#@<< Run parameters >>
#@+node:gcross.20090908092630.1905:<< Run parameters >>
for number_of_particles in [64]: #256]: #64,128,256]:
    for harmonic_oscillator_strength in [1]: #,10,100]:
        for hard_sphere_radius in [0.1]: #[0.01,0.1]:
            for backflow_coefficient in [1]: #0,0.01,0.1,1,-0.01,-0.1,-1]:
#@-node:gcross.20090908092630.1905:<< Run parameters >>
#@nl
                try:
                    #@                    << Run simulation for given parameters >>
                    #@+node:gcross.20090908092630.1864:<< Run simulation for given parameters >>
                    my_directory = "{output_root_directory}/{number_of_particles}/{harmonic_oscillator_strength}/{hard_sphere_radius}".format(**vars()) 
                    if (my_rank == 0):
                        print
                        print system_properties_message.format(**vars())

                    #@<< Final configuration >>
                    #@+node:gcross.20090908092630.1865:<< Final configuration >>
                    for parameter_name in [
                            "hard_sphere_radius",
                            "number_of_particles",
                            "backflow_coefficient",
                        ]:
                        configuration[parameter_name] = vars()[parameter_name]
                    configuration["number_of_rotating_particles"] = configuration["number_of_particles"]
                    configuration["harmonic_oscillator_coefficients"] = \
                        array([harmonic_oscillator_strength,]*configuration["number_of_dimensions"],dtype=double)
                    #@-node:gcross.20090908092630.1865:<< Final configuration >>
                    #@nl

                    system = System(**configuration)
                    system.rotation_rate = system.number_of_rotating_particles / system.number_of_particles
                    #@<< Initialize physics >>
                    #@+node:gcross.20090908092630.1866:<< Initialize physics >>
                    for physics in [
                        HarmonicOscillator,
                        RotationEffectivePotential,
                        SecondOrderGreensFunction,
                        HardSphereInteraction
                        ]: system.add_physics(physics)
                    #@-node:gcross.20090908092630.1866:<< Initialize physics >>
                    #@nl
                    #@<< Initialize observables >>
                    #@+node:gcross.20090908092630.1867:<< Initialize observables >>
                    for slice_name, slice_number in [("left",0),("center",system.center_slice_number)]: #,("right",system.number_of_slices-1)]:
                        density_slice_subdirectory = "{my_directory}/{backflow_coefficient}/{slice_name}".format(**vars())
                        for observable in [
                                #@            << Slice observables >>
                                #@+node:gcross.20090908092630.1868:<< Slice observables >>
                                RadialDensityHistogram(
                                    slice_number,
                                    radial_densities_histogram_maximum_radius,
                                    radial_densities_histogram_bin_count,
                                    density_slice_subdirectory + "/radial-density"
                                ),
                                PlaneRadialDensityHistogram(
                                    slice_number,
                                    radial_densities_histogram_maximum_radius,
                                    radial_densities_histogram_bin_count,
                                    density_slice_subdirectory + "/plane-radial-density"
                                ),
                                ParticleSeparationHistogram(
                                    slice_number,
                                    5,100,
                                    density_slice_subdirectory + "/particle-separation".format(**vars())
                                ),
                                #@-node:gcross.20090908092630.1868:<< Slice observables >>
                                #@nl
                            ]:
                            pass
                            system.add_observable(observable)

                    center_slice = system.number_of_slices // 2
                    for observable in  [
                        #@    << Averaged observables >>
                        #@+node:gcross.20090908092630.1908:<< Averaged observables >>
                        AverageParticleSeparationEstimate(
                            center_slice,
                            "{my_directory}/particle-separation".format(**vars()),
                            backflow_coefficient
                        ),
                        TotalEnergyEstimate(
                            "{my_directory}/total-energy".format(**vars()),
                            backflow_coefficient
                        ),
                        #@nonl
                        #@-node:gcross.20090908092630.1908:<< Averaged observables >>
                        #@nl
                        ]: system.add_observable(observable)
                    #@-node:gcross.20090908092630.1867:<< Initialize observables >>
                    #@nl
                    system.run()
                    system.total_and_write_observables()
                    del system.observables
                    del system
                    gc.collect()
                    #@-node:gcross.20090908092630.1864:<< Run simulation for given parameters >>
                    #@nl
                except FailureToContructValidSystemException, e:
                    print e
                    print "Aborting..."
                    comm.Abort(-1)
#@-node:gcross.20090908092630.1851:@thin run.py
#@-leo
