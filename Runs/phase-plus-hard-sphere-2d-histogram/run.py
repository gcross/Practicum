#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20100113111641.2614:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20100113111641.2616:<< Imports >>
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

from scipy.misc import derivative
import __builtin__
from itertools import imap, combinations
#@-node:gcross.20100113111641.2616:<< Imports >>
#@nl

#@+others
#@+node:gcross.20100113111641.2617:class RotationEffectivePotential
class RotationEffectivePotential(Physics):
    #@    @+others
    #@+node:gcross.20100113111641.2618:hooks
    hooks = ["effective_potentials"]
    #@-node:gcross.20100113111641.2618:hooks
    #@+node:gcross.20100113111641.2619:__init__
    def __init__(self,system):
        Physics.__init__(self,system)
        self.rotation_plane_axis_1 = system.rotation_plane_axis_1
        self.rotation_plane_axis_2 = system.rotation_plane_axis_2

        self.number_of_rotating_particles = system.number_of_rotating_particles
        self.frame_angular_velocity = system.frame_angular_velocity

        self.lambda_ = system.lambda_
    #@-node:gcross.20100113111641.2619:__init__
    #@+node:gcross.20100113111641.2620:accumulate_potential
    def accumulate_potential(self,x,xij2,U,gradU2):
        gradient_phase = vpif.angular_momentum.compute_gradient_ho_phase_choice(
            x,
            self.number_of_rotating_particles
        )

        vpif.angular_momentum.accumulate_rotating_frame_potential(
            x, gradient_phase,
            self.frame_angular_velocity,self.lambda_,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            U
        )
    #@-node:gcross.20100113111641.2620:accumulate_potential
    #@-others
#@-node:gcross.20100113111641.2617:class RotationEffectivePotential
#@+node:gcross.20100113111641.2621:class HarmonicOscillator
class HarmonicOscillator(Physics):
  #@  @+others
  #@+node:gcross.20100113111641.2622:hooks
  hooks = ["physical_potentials","trial_functions"]
  #@-node:gcross.20100113111641.2622:hooks
  #@+node:gcross.20100113111641.2623:__init__
  def __init__(self,system):
      Physics.__init__(self,system)
      try:
          self.potential_coefficients = system.potential_harmonic_oscillator_coefficients
          self.trial_coefficients = sqrt(system.trial_harmonic_oscillator_coefficients)
      except AttributeError:
              raise ValueError("System needs to define 'harmonic_oscillator_coefficients' to use harmonic oscillator physics!")       
  #@-node:gcross.20100113111641.2623:__init__
  #@+node:gcross.20100113111641.2624:accumulate_potential
  def accumulate_potential(self,x,xij2,U,gradU2):
      vpif.harmonic_oscillator.accumulate_potential(
          x,
          self.potential_coefficients,
          U,
          gradU2
      )
  #@-node:gcross.20100113111641.2624:accumulate_potential
  #@+node:gcross.20100113111641.2625:compute_trial_weight
  def compute_trial_weight(self,x,xij2):
      return vpif.harmonic_oscillator.compute_trial_weight(
          x,
          self.trial_coefficients
      )
  #@-node:gcross.20100113111641.2625:compute_trial_weight
  #@+node:gcross.20100113111641.2626:accumulate_trial_derivatives
  def accumulate_trial_derivatives(self,
          x,xij2,
          gradient_of_log_trial_fn,laplacian_of_log_trial_fn
      ):
      vpif.harmonic_oscillator.accumulate_trial_derivatives(
          x,self.trial_coefficients,
          gradient_of_log_trial_fn,laplacian_of_log_trial_fn
      )
  #@-node:gcross.20100113111641.2626:accumulate_trial_derivatives
  #@-others
#@-node:gcross.20100113111641.2621:class HarmonicOscillator
#@-others

#@<< System Configuration >>
#@+node:gcross.20100113111641.2615:<< System Configuration >>
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
    "total_number_of_observations": 1000,
    "number_of_prethermalization_steps": 1000,
    # Move parameters
    "dM": 80,
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.01,0.5,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
}
configuration["trial_harmonic_oscillator_coefficients"] = \
    array([1.0,]*configuration["number_of_dimensions"],dtype=double)
configuration["potential_harmonic_oscillator_coefficients"] = \
    array([1.0,]*configuration["number_of_dimensions"],dtype=double)
#@-node:gcross.20100113111641.2615:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20100113111641.2627:<< Histogram Configuration >>
_1d_densities_histogram_left  = [-2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_right = [+2,]*configuration["number_of_dimensions"]
_1d_densities_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 50
#@-node:gcross.20100113111641.2627:<< Histogram Configuration >>
#@nl
#@<< System properties message >>
#@+node:gcross.20100113111641.2628:<< System properties message >>
system_properties_message = """\
Examining system with
    *) a frame rotating with velocity {frame_angular_velocity};
    *) a hard sphere radius of {hard_sphere_radius};
    *) {number_of_particles} particles;
    *) and {number_of_rotating_particles} particles rotating"""
#@-node:gcross.20100113111641.2628:<< System properties message >>
#@nl

output_root_directory = sys.argv[1]

for frame_angular_velocity in [0.7]:
    for hard_sphere_radius in [0]: #,0.1,0.2]:
        for number_of_particles in [6]:
            for number_of_rotating_particles in xrange(number_of_particles+1):
                #@                << Run simulation for given parameters >>
                #@+node:gcross.20100113111641.2629:<< Run simulation for given parameters >>
                my_directory = "{output_root_directory}/{frame_angular_velocity}/{hard_sphere_radius}/{number_of_particles}".format(**vars()) 
                if (my_rank == 0):
                    print
                    print system_properties_message.format(**vars())

                for parameter_name in ["hard_sphere_radius","number_of_particles","frame_angular_velocity","number_of_rotating_particles"]:
                    configuration[parameter_name] = vars()[parameter_name]

                system = System(**configuration)
                #@<< Initialize physics >>
                #@+node:gcross.20100113111641.2630:<< Initialize physics >>
                for physics in [
                    HarmonicOscillator,
                    RotationEffectivePotential,
                    SecondOrderGreensFunction,
                    HardSphereInteraction,
                    ]: system.add_physics(physics)
                #@-node:gcross.20100113111641.2630:<< Initialize physics >>
                #@nl
                #@<< Initialize observables >>
                #@+node:gcross.20100113111641.2631:<< Initialize observables >>
                for slice_name, slice_number in [("left",0),("center",system.center_slice_number),("right",system.number_of_slices-1)]:
                    density_slice_subdirectory = "{my_directory}/{number_of_rotating_particles}/{slice_name}".format(**vars())
                    for observable in [
                            PositionDensity2DHistogram(
                                slice_number,
                                -2,-2,
                                +2,+2,
                                20,
                                "{density_slice_subdirectory}/density-2d".format(**vars())
                            ),
                        ]:
                        system.add_observable(observable)

                center_slice = system.number_of_slices // 2
                for observable in  [
                    TotalEnergyEstimate("{my_directory}/total-energy".format(**vars()),number_of_rotating_particles),
                    ]: system.add_observable(observable)
                #@-node:gcross.20100113111641.2631:<< Initialize observables >>
                #@nl
                system.run()
                system.total_and_write_observables()
                del system.observables
                del system
                gc.collect()
                #@-node:gcross.20100113111641.2629:<< Run simulation for given parameters >>
                #@nl
#@-node:gcross.20100113111641.2614:@thin run.py
#@-leo
