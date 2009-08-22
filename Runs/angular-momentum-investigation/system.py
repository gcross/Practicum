#@+leo-ver=4-thin
#@+node:gcross.20090821120721.1291:@thin system.py
import sys
sys.path.append("../../lib")

#@<< Imports >>
#@+node:gcross.20090821120721.1292:<< Imports >>
import vpi

from numpy import *
from numpy.random import rand

import os
import os.path

from itertools import izip
#@-node:gcross.20090821120721.1292:<< Imports >>
#@nl

#@<< MPI Initialization >>
#@+node:gcross.20090821120721.1294:<< MPI Initialization >>
from mpi4py import MPI

comm = MPI.COMM_WORLD

number_of_processors = comm.Get_size()
my_rank = comm.Get_rank()
#@-node:gcross.20090821120721.1294:<< MPI Initialization >>
#@nl

#@+others
#@+node:gcross.20090821120721.1310:Observable classes
#@+others
#@+node:gcross.20090821174249.1325:class Observable
class Observable(object):
    #@    @+others
    #@+node:gcross.20090821174249.1326:total_and_write
    def total_and_write(self):
        totals = self.compute_total()
        if my_rank == 0:
            self.write_out_totals(totals)
    #@-node:gcross.20090821174249.1326:total_and_write
    #@-others
#@-node:gcross.20090821174249.1325:class Observable
#@+node:gcross.20090821120721.1325:Histograms
#@+node:gcross.20090821120721.1316:class Histogram
class Histogram(Observable):
    #@    @+others
    #@+node:gcross.20090821120721.1317:compute_total
    def compute_total(self):
        total_histogram = zeros(self.histogram.shape,dtype='i',order='Fortran')
        comm.Reduce((self.histogram,MPI.INT),(total_histogram,MPI.INT))
        return total_histogram
    #@-node:gcross.20090821120721.1317:compute_total
    #@-others
#@-node:gcross.20090821120721.1316:class Histogram
#@+node:gcross.20090821120721.1311:class PositionDensity1DHistogram
class PositionDensity1DHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090821120721.1312:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filenames):
        assert len(left) == len(right)
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((len(left),number_of_bins),dtype='i',order='Fortran')
        self.filenames = filenames
    #@-node:gcross.20090821120721.1312:__init__
    #@+node:gcross.20090821120721.1313:update
    def update(self):
        vpi.histograms.accumulate_1d_densities(
            self.system.x[self.slice_number],
            self.left,self.right,
            self.histogram
        )
    #@-node:gcross.20090821120721.1313:update
    #@+node:gcross.20090821120721.1326:write_out_totals
    def write_out_totals(self,histograms):
        for filename, histogram, left, right in izip(self.filenames,histograms,self.left,self.right):
            ensure_path_to_file_exists(filename)
            with open(filename,"w") as f:
                total_counts = float(sum(histogram))
                bin_width = float(right-left)/self.number_of_bins
                current = float(left)+bin_width/2
                for count in histogram:
                    print >> f, "{0} {1}".format(current,count/total_counts)
                    current += bin_width
    #@-node:gcross.20090821120721.1326:write_out_totals
    #@-others
#@-node:gcross.20090821120721.1311:class PositionDensity1DHistogram
#@+node:gcross.20090821120721.1322:class RadialDensityHistogram
class RadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090821120721.1323:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090821120721.1323:__init__
    #@+node:gcross.20090821120721.1324:update
    def update(self):
        vpi.histograms.accumulate_radial_densities(
            self.system.x[self.slice_number],
            self.maximum_radius,
            self.histogram
        )
    #@-node:gcross.20090821120721.1324:update
    #@+node:gcross.20090821120721.1328:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(self.maximum_radius)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                normalization = total_counts*4.0/3.0*pi*((current+bin_width/2)**3-(current-bin_width/2)**3)
                print >> f, "{0} {1}".format(current,count/normalization)
                current += bin_width
    #@-node:gcross.20090821120721.1328:write_out_totals
    #@-others
#@-node:gcross.20090821120721.1322:class RadialDensityHistogram
#@-node:gcross.20090821120721.1325:Histograms
#@+node:gcross.20090821174249.1328:Energy estimates
#@+node:gcross.20090821174249.1330:class TotalEnergyEstimate
class TotalEnergyEstimate(Observable):
    #@    @+others
    #@+node:gcross.20090821174249.1331:__init__
    def __init__(self,filename,label):
        self.estimates = [0,0]
        self.filename = filename
        self.label = label
    #@-node:gcross.20090821174249.1331:__init__
    #@+node:gcross.20090821174249.1332:update
    def update(self):
        system = self.system
        for i, slice_number in enumerate([0,-1]):
            gradient_of_log_trial_fn, laplacian_of_log_trial_fn = system.compute_trial_derivatives(system.x[i],system.xij2[i])
            self.estimates[i] += \
                vpi.observables.compute_local_energy_estimate(
                    system.U[slice_number],
                    gradient_of_log_trial_fn, laplacian_of_log_trial_fn,
                    system.lambda_,
                )
    #@-node:gcross.20090821174249.1332:update
    #@+node:gcross.20090821174249.1338:compute_total
    def compute_total(self):
        total_estimates = zeros((len(self.estimates),),dtype=double)
        comm.Reduce((array(self.estimates,dtype=double),MPI.DOUBLE),(total_estimates,MPI.DOUBLE))
        total_estimates /= self.system.total_number_of_observations
        return total_estimates
    #@-node:gcross.20090821174249.1338:compute_total
    #@+node:gcross.20090821174249.1336:write_out_totals
    def write_out_totals(self,total_estimates):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"a") as f:
            print >> f, self.label, mean(total_estimates), std(total_estimates)
    #@-node:gcross.20090821174249.1336:write_out_totals
    #@-others
#@-node:gcross.20090821174249.1330:class TotalEnergyEstimate
#@+node:gcross.20090821174249.1335:class SliceEnergyEstimate
class SliceEnergyEstimate(Observable):
    #@    @+others
    #@+node:gcross.20090821174249.1339:(fields)
    estimate = 0
    #@-node:gcross.20090821174249.1339:(fields)
    #@+node:gcross.20090821174249.1340:__init__
    def __init__(self,slice_number,filename,label):
        self.slice_number = slice_number
        self.filename = filename
        self.label = label
    #@-node:gcross.20090821174249.1340:__init__
    #@+node:gcross.20090821174249.1343:compute_total
    def compute_total(self):
        estimate = zeros((1,),dtype='d',order='Fortran')
        estimate[0] = self.estimate
        total_estimate = zeros((1,),dtype='d',order='Fortran')
        comm.Reduce((estimate,MPI.DOUBLE),(total_estimate,MPI.DOUBLE))
        return total_estimate[0] / self.system.total_number_of_observations
    #@-node:gcross.20090821174249.1343:compute_total
    #@+node:gcross.20090821174249.1345:write_out_totals
    def write_out_totals(self,total_estimate):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"a") as f:
            print >> f, self.label, total_estimate
    #@-node:gcross.20090821174249.1345:write_out_totals
    #@-others
#@-node:gcross.20090821174249.1335:class SliceEnergyEstimate
#@+node:gcross.20090821174249.1346:class EffectivePotentialEnergyEstimate
class EffectivePotentialEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090821174249.1347:update
    def update(self):
        system = self.system
        U = zeros((1,system.number_of_particles),dtype=double,order='Fortran')
        gradU = zeros((1,system.number_of_particles,system.number_of_dimensions),dtype=double,order='Fortran')
        vpi.angular_momentum.compute_effective_rotational_potential(
            system.x[self.slice_number:self.slice_number+1],system.lambda_,
            system.rotation_plane_axis_1,system.rotation_plane_axis_2,
            system.frame_angular_velocity,system.number_of_rotating_particles,
            U, gradU
        )
        self.estimate += sum(U)
    #@-node:gcross.20090821174249.1347:update
    #@-others
#@-node:gcross.20090821174249.1346:class EffectivePotentialEnergyEstimate
#@+node:gcross.20090821174249.1350:class PhysicalPotentialEnergyEstimate
class PhysicalPotentialEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090821174249.1351:update
    def update(self):
        self.estimate += sum(dot(self.system.x[self.slice_number]**2,self.system.harmonic_oscillator_coefficients))/2.0
    #@-node:gcross.20090821174249.1351:update
    #@-others
#@-node:gcross.20090821174249.1350:class PhysicalPotentialEnergyEstimate
#@+node:gcross.20090821174249.1352:class TotalPotentialEnergyEstimate
class TotalPotentialEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090821174249.1354:update
    def update(self):
        self.estimate += sum(self.system.U[self.slice_number])
    #@-node:gcross.20090821174249.1354:update
    #@-others
#@-node:gcross.20090821174249.1352:class TotalPotentialEnergyEstimate
#@-node:gcross.20090821174249.1328:Energy estimates
#@-others
#@-node:gcross.20090821120721.1310:Observable classes
#@+node:gcross.20090821164402.1307:Functions
#@+node:gcross.20090821164402.1308:ensure_path_to_file_exists
def ensure_path_to_file_exists(path):
    directory, _ = os.path.split(path)
    if not os.path.exists(directory):
        os.makedirs(directory)
#@-node:gcross.20090821164402.1308:ensure_path_to_file_exists
#@-node:gcross.20090821164402.1307:Functions
#@+node:gcross.20090821174249.1315:class System
class System:
    #@    @+others
    #@+node:gcross.20090821174249.1318:Physics Functions
    #@+node:gcross.20090821120721.1298:compute_potential
    def compute_potential(self,x,_):
        x_sq = x**2
        U = array(dot(x_sq,self.harmonic_oscillator_coefficients)/2.0,dtype=double,order='Fortran')
        gradU  = array(x,dtype=double,order='Fortran')
        gradU *= self.harmonic_oscillator_coefficients
        vpi.angular_momentum.compute_effective_rotational_potential(
            x,self.lambda_,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            self.frame_angular_velocity,self.number_of_rotating_particles,
            U, gradU
        )
        gradU2 = sum(sum(gradU**2,axis=-1),axis=-1)
        return U, gradU2, False
    #@-node:gcross.20090821120721.1298:compute_potential
    #@+node:gcross.20090821120721.1299:compute_trial
    def compute_trial(self,x,_):
        return -sum(dot(x**2,self.harmonic_oscillator_coefficients))/2
    #@-node:gcross.20090821120721.1299:compute_trial
    #@+node:gcross.20090821120721.1300:compute_trial_derivatives
    def compute_trial_derivatives(self,x,xij2):
        gradient_of_log_trial_fn = -x*self.harmonic_oscillator_coefficients
        laplacian_of_log_trial_fn = -sum(self.harmonic_oscillator_coefficients)*x.shape[0]
        return gradient_of_log_trial_fn, laplacian_of_log_trial_fn
    #@-node:gcross.20090821120721.1300:compute_trial_derivatives
    #@-node:gcross.20090821174249.1318:Physics Functions
    #@+node:gcross.20090821174249.1322:Observable management
    #@+node:gcross.20090821174249.1323:add_observable
    def add_observable(self,observable):
        self.observables.append(observable)
        observable.system = self
    #@-node:gcross.20090821174249.1323:add_observable
    #@+node:gcross.20090821174249.1324:total_and_write_observables
    def total_and_write_observables(self):
        for observable in self.observables:
            observable.total_and_write()
    #@-node:gcross.20090821174249.1324:total_and_write_observables
    #@-node:gcross.20090821174249.1322:Observable management
    #@+node:gcross.20090821174249.1316:__init__
    def __init__(self,**keywords):
        self.__dict__.update(keywords)

        number_of_slices = self.number_of_slices
        number_of_particles = self.number_of_particles
        number_of_dimensions = self.number_of_dimensions

        vpi.rand_utils.init_seed(my_rank)

        self.x = vpi.lattice.make_lattice(1.0,number_of_slices,number_of_particles,number_of_dimensions)

        self.xij2 = zeros((number_of_slices,number_of_particles,number_of_particles),dtype=double,order='Fortran')
        vpi.xij.update_xij(self.xij2,self.x)

        self.U = zeros((number_of_slices,number_of_particles),dtype=double,order='Fortran')
        self.gradU2 = zeros((number_of_slices),dtype=double,order='Fortran')

        self.slice_move_attempted_counts = zeros((number_of_slices,),'i')
        self.slice_move_accepted_counts = zeros((number_of_slices,),'i')

        self.move_type_attempted_counts = zeros((3,),'i')
        self.move_type_accepted_counts = zeros((3,),'i')

        U_weights, gU2_weights = vpi.gfn.initialize_4th_order_weights(number_of_slices)
        self.U_weights = U_weights
        self.gU2_weights = gU2_weights

        self.number_of_observations = self.total_number_of_observations // number_of_processors + 1
        self.total_number_of_observations = self.number_of_observations * number_of_processors

        self.number_of_thermalizations_per_observation = number_of_particles * number_of_slices // self.dM

        assert (number_of_slices % 2 == 0 and number_of_slices % 4 == 2)
        self.center_slice_number = number_of_slices // 2

        self.observables = []

    #@-node:gcross.20090821174249.1316:__init__
    #@+node:gcross.20090821120721.1304:run
    def run(self):
        #@    << Stash properties into local variables >>
        #@+node:gcross.20090821174249.1319:<< Stash properties into local variables >>
        x = self.x
        xij2 = self.xij2
        U = self.U
        gradU2 = self.gradU2
        number_of_prethermalization_steps = self.number_of_prethermalization_steps
        number_of_thermalizations_per_observation = self.number_of_thermalizations_per_observation
        move_type_probabilities = self.move_type_probabilities
        move_type_differentials = self.move_type_differentials
        dM = self.dM
        lambda_ = self.lambda_
        low_swap_dimension = self.low_swap_dimension
        high_swap_dimension = self.high_swap_dimension
        slice_move_attempted_counts = self.slice_move_attempted_counts
        move_type_attempted_counts = self.move_type_attempted_counts
        slice_move_accepted_counts = self.slice_move_accepted_counts
        move_type_accepted_counts = self.move_type_accepted_counts
        compute_potential = self.compute_potential
        compute_trial = self.compute_trial
        U_weights = self.U_weights
        gU2_weights = self.gU2_weights
        use_4th_order_green_function = self.use_4th_order_green_function

        observables = self.observables
        #@-node:gcross.20090821174249.1319:<< Stash properties into local variables >>
        #@nl
        #@    << Prethermalize the system >>
        #@+node:gcross.20090821174249.1320:<< Prethermalize the system >>
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
            compute_potential,compute_trial,
            U_weights,gU2_weights,
            use_4th_order_green_function,
        )
        #@nonl
        #@-node:gcross.20090821174249.1320:<< Prethermalize the system >>
        #@nl
        #@    << Main iteration >>
        #@+node:gcross.20090821174249.1321:<< Main iteration >>
        decile = self.number_of_observations // 10
        for number_completed in xrange(self.number_of_observations):
            vpi.thermalize.thermalize_path(
                x,xij2,
                U,gradU2,
                number_of_thermalizations_per_observation,
                move_type_probabilities,move_type_differentials,
                dM,
                lambda_,
                low_swap_dimension,high_swap_dimension,
                slice_move_attempted_counts,move_type_attempted_counts,
                slice_move_accepted_counts,move_type_accepted_counts,
                compute_potential,compute_trial,
                U_weights,gU2_weights,
                use_4th_order_green_function,
            )
            for observable in observables:
                observable.update()
            if (number_completed % decile == 0) and (my_rank == 0):
                print "{0:.0%} complete;  local bridge move acceptance rate = {1:.0%}, local rigid move acceptance rate = {2:.0%}".format(
                    float(number_completed)/self.number_of_observations,
                    float(move_type_accepted_counts[0])/move_type_attempted_counts[0],
                    float(move_type_accepted_counts[1])/move_type_attempted_counts[1],
                )
        #@-node:gcross.20090821174249.1321:<< Main iteration >>
        #@nl
    #@-node:gcross.20090821120721.1304:run
    #@-others
#@-node:gcross.20090821174249.1315:class System
#@-others
#@-node:gcross.20090821120721.1291:@thin system.py
#@-leo
