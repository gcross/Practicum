#@+leo-ver=4-thin
#@+node:gcross.20090827130017.1614:@thin system.py
import sys
sys.path.append("../../lib")

#@<< Imports >>
#@+node:gcross.20090827130017.1615:<< Imports >>
import vpi

from numpy import *
from numpy.random import rand

import os
import os.path

import itertools
from itertools import izip
#@-node:gcross.20090827130017.1615:<< Imports >>
#@nl

#@<< MPI Initialization >>
#@+node:gcross.20090827130017.1616:<< MPI Initialization >>
from mpi4py import MPI

comm = MPI.COMM_WORLD

number_of_processors = comm.Get_size()
my_rank = comm.Get_rank()
#@-node:gcross.20090827130017.1616:<< MPI Initialization >>
#@nl

#@+others
#@+node:gcross.20090827130017.1617:Observable classes
#@+others
#@+node:gcross.20090827130017.1762:Base classes
#@+node:gcross.20090827130017.1618:class Observable
class Observable(object):
    #@    @+others
    #@+node:gcross.20090827130017.1619:total_and_write
    def total_and_write(self):
        totals = self.compute_total()
        if my_rank == 0:
            self.write_out_totals(totals)
    #@-node:gcross.20090827130017.1619:total_and_write
    #@-others
#@-node:gcross.20090827130017.1618:class Observable
#@+node:gcross.20090827130017.1620:class AverageValuesEstimate
class AverageValuesEstimate(Observable):
    #@    @+others
    #@+node:gcross.20090827130017.1621:compute_total
    def compute_total(self):
        total_estimate = zeros(self.estimate.shape,dtype='d',order='Fortran')
        comm.Reduce((self.estimate,MPI.DOUBLE),(total_estimate,MPI.DOUBLE))
        total_estimate /= self.system.total_number_of_observations
        return total_estimate
    #@-node:gcross.20090827130017.1621:compute_total
    #@-others
#@-node:gcross.20090827130017.1620:class AverageValuesEstimate
#@+node:gcross.20090827130017.1622:class SingleAverageValueEstimate
class SingleAverageValueEstimate(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1623:__init__
    def __init__(self):
        self.estimate = zeros((1,),dtype=double)
    #@-node:gcross.20090827130017.1623:__init__
    #@-others
#@-node:gcross.20090827130017.1622:class SingleAverageValueEstimate
#@+node:gcross.20090827130017.1624:class EstimatesAppendedToFile
class EstimatesAppendedToFile(Observable):
    #@    @+others
    #@+node:gcross.20090827130017.1625:__init__
    def __init__(self,filename,label):
        self.filename = filename
        self.label = label
    #@-node:gcross.20090827130017.1625:__init__
    #@+node:gcross.20090827130017.1626:write_out_totals
    def write_out_totals(self,total_estimate):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"a") as f:
            print >> f, self.label,
            for value in total_estimate:
                print >> f, value,
            print >> f

    #@-node:gcross.20090827130017.1626:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1624:class EstimatesAppendedToFile
#@+node:gcross.20090827130017.1627:class SingleAverageValueEstimateAppendedToFile
class SingleAverageValueEstimateAppendedToFile(SingleAverageValueEstimate,EstimatesAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1628:__init__
    def __init__(self,filename,label):
        SingleAverageValueEstimate.__init__(self)
        EstimatesAppendedToFile.__init__(self,filename,label)
    #@-node:gcross.20090827130017.1628:__init__
    #@-others
#@-node:gcross.20090827130017.1627:class SingleAverageValueEstimateAppendedToFile
#@+node:gcross.20090827130017.1629:class SingleAverageValueAtSliceEstimateAppendedToFile
class SingleAverageValueAtSliceEstimateAppendedToFile(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1630:__init__
    def __init__(self,slice_number,label,filename):
        SingleAverageValueEstimateAppendedToFile.__init__(self,label,filename)
        self.slice_number = slice_number
    #@-node:gcross.20090827130017.1630:__init__
    #@-others
#@-node:gcross.20090827130017.1629:class SingleAverageValueAtSliceEstimateAppendedToFile
#@-node:gcross.20090827130017.1762:Base classes
#@+node:gcross.20090827130017.1631:Histograms
#@+node:gcross.20090827130017.1632:class Histogram
class Histogram(Observable):
    #@    @+others
    #@+node:gcross.20090827130017.1633:compute_total
    def compute_total(self):
        total_histogram = zeros(self.histogram.shape,dtype='i',order='Fortran')
        comm.Reduce((self.histogram,MPI.INT),(total_histogram,MPI.INT))
        return total_histogram
    #@-node:gcross.20090827130017.1633:compute_total
    #@+node:gcross.20090827130017.1634:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram))
        bin_width = float(self.right-self.left)/self.number_of_bins
        current = float(self.left)+bin_width/2
        with open(self.filename,"w") as f:
            for count in histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090827130017.1634:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1632:class Histogram
#@+node:gcross.20090827130017.1635:class PositionDensity1DHistogram
class PositionDensity1DHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1636:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filenames):
        assert len(left) == len(right)
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((len(left),number_of_bins),dtype='i',order='Fortran')
        self.filenames = filenames
    #@-node:gcross.20090827130017.1636:__init__
    #@+node:gcross.20090827130017.1637:update
    def update(self):
        vpi.histograms.accumulate_1d_densities(
            self.system.x[self.slice_number],
            self.left,self.right,
            self.histogram
        )
    #@-node:gcross.20090827130017.1637:update
    #@+node:gcross.20090827130017.1638:write_out_totals
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
    #@-node:gcross.20090827130017.1638:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1635:class PositionDensity1DHistogram
#@+node:gcross.20090827130017.1639:class RadialDensityHistogram
class RadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1640:(fields)
    left = 0
    #@-node:gcross.20090827130017.1640:(fields)
    #@+node:gcross.20090827130017.1641:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1641:__init__
    #@+node:gcross.20090827130017.1642:update
    def update(self):
        vpi.histograms.accumulate_radial_densities(
            self.system.x[self.slice_number],
            self.maximum_radius,
            self.histogram
        )
    #@-node:gcross.20090827130017.1642:update
    #@-others
#@-node:gcross.20090827130017.1639:class RadialDensityHistogram
#@+node:gcross.20090827130017.1643:class PlaneRadialDensityHistogram
class PlaneRadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1644:(fields)
    left = 0
    #@-node:gcross.20090827130017.1644:(fields)
    #@+node:gcross.20090827130017.1645:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1645:__init__
    #@+node:gcross.20090827130017.1646:update
    def update(self):
        vpi.histograms.accumulate_plane_radial_densities(
            self.system.x[self.slice_number],
            self.maximum_radius,
            self.system.rotation_plane_axis_1,
            self.system.rotation_plane_axis_2,
            self.histogram
        )
    #@-node:gcross.20090827130017.1646:update
    #@-others
#@-node:gcross.20090827130017.1643:class PlaneRadialDensityHistogram
#@+node:gcross.20090827130017.1647:class RecipricalRadiusSquaredDensityHistogram
class RecipricalRadiusSquaredDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1648:(fields)
    left = 0
    #@-node:gcross.20090827130017.1648:(fields)
    #@+node:gcross.20090827130017.1649:__init__
    def __init__(self,slice_number,maximum_value,number_of_bins,filename):
        self.right = self.maximum_value = maximum_value
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1649:__init__
    #@+node:gcross.20090827130017.1650:update
    def update(self):
        vpi.histograms.accumulate_reciprical_radius_squared_densities(
            self.system.x[self.slice_number],
            self.maximum_value,
            self.histogram
        )
    #@-node:gcross.20090827130017.1650:update
    #@-others
#@-node:gcross.20090827130017.1647:class RecipricalRadiusSquaredDensityHistogram
#@+node:gcross.20090827130017.1651:class RecipricalPlaneRadiusSquaredDensityHistogram
class RecipricalPlaneRadiusSquaredDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1652:(fields)
    left = 0
    #@-node:gcross.20090827130017.1652:(fields)
    #@+node:gcross.20090827130017.1653:__init__
    def __init__(self,slice_number,maximum_value,number_of_bins,filename):
        self.right = self.maximum_value = maximum_value
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1653:__init__
    #@+node:gcross.20090827130017.1654:update
    def update(self):
        vpi.histograms.accumulate_recip_plane_r_sq_densities(
            self.system.x[self.slice_number],
            self.maximum_value,
            self.system.rotation_plane_axis_1,
            self.system.rotation_plane_axis_2,
            self.histogram
        )
    #@-node:gcross.20090827130017.1654:update
    #@+node:gcross.20090827130017.1655:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(self.maximum_value)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090827130017.1655:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1651:class RecipricalPlaneRadiusSquaredDensityHistogram
#@+node:gcross.20090827130017.1656:class AngularSeparationDensityHistogram
class AngularSeparationDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1657:(fields)
    left = 0
    right = 2*pi
    #@-node:gcross.20090827130017.1657:(fields)
    #@+node:gcross.20090827130017.1658:__init__
    def __init__(self,slice_number,number_of_bins,filename):
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1658:__init__
    #@+node:gcross.20090827130017.1659:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        vpi.histograms.accumulate_angular_separation_densities(
            angles,
            self.histogram
        )
    #@-node:gcross.20090827130017.1659:update
    #@+node:gcross.20090827130017.1660:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(2*pi)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090827130017.1660:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1656:class AngularSeparationDensityHistogram
#@+node:gcross.20090827130017.1661:class NeighborAngularSeparationDensityHistogram
class NeighborAngularSeparationDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1662:(fields)
    left = 0
    right = 2*pi
    #@-node:gcross.20090827130017.1662:(fields)
    #@+node:gcross.20090827130017.1663:__init__
    def __init__(self,slice_number,number_of_bins,filename):
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1663:__init__
    #@+node:gcross.20090827130017.1664:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        vpi.histograms.accumulate_neighbor_angular_separation_densities(
            angles,
            self.histogram
        )
    #@-node:gcross.20090827130017.1664:update
    #@-others
#@-node:gcross.20090827130017.1661:class NeighborAngularSeparationDensityHistogram
#@+node:gcross.20090827130017.1665:class AngularVelocityHistogram
class AngularVelocityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1666:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1666:__init__
    #@+node:gcross.20090827130017.1667:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        vpi.histograms.accumulate(first_derivatives,self.left,self.right,self.histogram)
    #@-node:gcross.20090827130017.1667:update
    #@-others
#@-node:gcross.20090827130017.1665:class AngularVelocityHistogram
#@+node:gcross.20090827130017.1668:class AngularVelocitySquaredHistogram
class AngularVelocitySquaredHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1669:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1669:__init__
    #@+node:gcross.20090827130017.1670:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        vpi.histograms.accumulate(first_derivatives,self.left,self.right,self.histogram)
    #@-node:gcross.20090827130017.1670:update
    #@-others
#@-node:gcross.20090827130017.1668:class AngularVelocitySquaredHistogram
#@+node:gcross.20090827130017.1671:class RotationQuadraticTermHistogram
class RotationQuadraticTermHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1672:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1672:__init__
    #@+node:gcross.20090827130017.1673:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            x,
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        term = (first_derivatives ** 2) / (x[:,system.rotation_plane_axis_2-1]**2+x[:,system.rotation_plane_axis_1-1]**2)
        vpi.histograms.accumulate(term,self.left,self.right,self.histogram)
    #@-node:gcross.20090827130017.1673:update
    #@-others
#@-node:gcross.20090827130017.1671:class RotationQuadraticTermHistogram
#@+node:gcross.20090827130017.1674:class AngularVelocityAndRadiusHistogram
class AngularVelocityAndRadiusHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1675:__init__
    def __init__(self,slice_number,maximum_angular_velocity,maximum_radius,number_of_bins,filename):
        self.maximum_angular_velocity = maximum_angular_velocity
        self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,)*2,dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1675:__init__
    #@+node:gcross.20090827130017.1676:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            x,
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        radii = sqrt(x[system.rotation_plane_axis_2-1]**2+x[system.rotation_plane_axis_1-1]**2)

        i_values = floor(first_derivatives/self.maximum_angular_velocity*self.number_of_bins)
        j_values = floor(radii/self.maximum_radius*self.number_of_bins)
        for (i,j) in izip(i_values,j_values):
            if (i >= 0) and (i < self.number_of_bins) and (j >= 0) and (j < self.number_of_bins):
                self.histogram[i,j] += 1
    #@-node:gcross.20090827130017.1676:update
    #@+node:gcross.20090827130017.1677:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram))
        with open(self.filename,"w") as f:
            for i in xrange(self.number_of_bins):
                for j in xrange(self.number_of_bins):
                    angular_velocity = (i+0.5)/self.number_of_bins*self.maximum_angular_velocity
                    radius = (j+0.5)/self.number_of_bins*self.maximum_radius
                    print >> f, "{0} {1} {2}".format(angular_velocity,radius,histogram[i,j]/total_counts)
                print >> f, ""
    #@-node:gcross.20090827130017.1677:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1674:class AngularVelocityAndRadiusHistogram
#@-node:gcross.20090827130017.1631:Histograms
#@+node:gcross.20090827130017.1678:Energy estimates
#@+node:gcross.20090827130017.1679:class TotalEnergyEstimate
class TotalEnergyEstimate(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1680:__init__
    def __init__(self,filename,label):
        self.filename = filename
        self.label = label
        self.estimate = zeros((2,),dtype=double)
    #@-node:gcross.20090827130017.1680:__init__
    #@+node:gcross.20090827130017.1681:update
    def update(self):
        system = self.system
        for i, slice_number in enumerate([0,-1]):
            gradient_of_log_trial_fn, laplacian_of_log_trial_fn = system.compute_trial_derivatives(system.x[i],system.xij2[i])
            self.estimate[i] += \
                vpi.observables.compute_local_energy_estimate(
                    system.U[slice_number],
                    gradient_of_log_trial_fn, laplacian_of_log_trial_fn,
                    system.lambda_,
                )
    #@-node:gcross.20090827130017.1681:update
    #@+node:gcross.20090827130017.1682:write_out_totals
    def write_out_totals(self,total_estimates):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"a") as f:
            print >> f, self.label, mean(total_estimates), std(total_estimates)
    #@-node:gcross.20090827130017.1682:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1679:class TotalEnergyEstimate
#@+node:gcross.20090827130017.1683:Slice estimates
#@+node:gcross.20090827130017.1684:class SliceEnergyEstimate
class SliceEnergyEstimate(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1685:__init__
    def __init__(self,slice_number,filename,label):
        SingleAverageValueEstimateAppendedToFile.__init__(self,filename,label)
        self.slice_number = slice_number
    #@-node:gcross.20090827130017.1685:__init__
    #@-others
#@-node:gcross.20090827130017.1684:class SliceEnergyEstimate
#@+node:gcross.20090827130017.1686:class EffectivePotentialSliceEnergyEstimate
class EffectivePotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1687:update
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
    #@-node:gcross.20090827130017.1687:update
    #@-others
#@-node:gcross.20090827130017.1686:class EffectivePotentialSliceEnergyEstimate
#@+node:gcross.20090827130017.1688:class PhysicalPotentialSliceEnergyEstimate
class PhysicalPotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1689:update
    def update(self):
        self.estimate += sum(dot(self.system.x[self.slice_number]**2,self.system.harmonic_oscillator_coefficients))/2.0
    #@-node:gcross.20090827130017.1689:update
    #@-others
#@-node:gcross.20090827130017.1688:class PhysicalPotentialSliceEnergyEstimate
#@+node:gcross.20090827130017.1690:class TotalPotentialSliceEnergyEstimate
class TotalPotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1691:update
    def update(self):
        self.estimate += sum(self.system.U[self.slice_number])
    #@-node:gcross.20090827130017.1691:update
    #@-others
#@-node:gcross.20090827130017.1690:class TotalPotentialSliceEnergyEstimate
#@-node:gcross.20090827130017.1683:Slice estimates
#@+node:gcross.20090827130017.1692:Path estimates
#@+node:gcross.20090827130017.1693:class PathEnergyEstimates
class PathEnergyEstimates(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1694:(fields)
    estimates = 0
    #@-node:gcross.20090827130017.1694:(fields)
    #@+node:gcross.20090827130017.1695:__init__
    def __init__(self,filename):
        self.filename = filename
    #@-node:gcross.20090827130017.1695:__init__
    #@+node:gcross.20090827130017.1696:write_out_totals
    def write_out_totals(self,total_estimates):
        ensure_path_to_file_exists(self.filename)
        center_slice_number = self.system.center_slice_number
        with open(self.filename,"w") as f:
            for slice_number, estimate in enumerate(total_estimates):
                print >> f, center_slice_number-abs(center_slice_number-slice_number), estimate, slice_number
    #@-node:gcross.20090827130017.1696:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1693:class PathEnergyEstimates
#@+node:gcross.20090827130017.1697:class EffectivePotentialPathEnergyEstimates
class EffectivePotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090827130017.1698:update
    def update(self):
        system = self.system
        U = zeros((system.number_of_slices,system.number_of_particles),dtype=double,order='Fortran')
        gradU = zeros((system.number_of_slices,system.number_of_particles,system.number_of_dimensions),dtype=double,order='Fortran')
        vpi.angular_momentum.compute_effective_rotational_potential(
            system.x,system.lambda_,
            system.rotation_plane_axis_1,system.rotation_plane_axis_2,
            system.frame_angular_velocity,system.number_of_rotating_particles,
            U, gradU
        )
        self.estimates += sum(U,axis=-1)
    #@-node:gcross.20090827130017.1698:update
    #@-others
#@-node:gcross.20090827130017.1697:class EffectivePotentialPathEnergyEstimates
#@+node:gcross.20090827130017.1699:class PhysicalPotentialPathEnergyEstimates
class PhysicalPotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090827130017.1700:update
    def update(self):
        self.estimates += sum(dot(self.system.x**2,self.system.harmonic_oscillator_coefficients),axis=-1)/2.0
    #@-node:gcross.20090827130017.1700:update
    #@-others
#@-node:gcross.20090827130017.1699:class PhysicalPotentialPathEnergyEstimates
#@+node:gcross.20090827130017.1701:class TotalPotentialPathEnergyEstimates
class TotalPotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090827130017.1702:update
    def update(self):
        self.estimates += sum(self.system.U,axis=-1)
    #@-node:gcross.20090827130017.1702:update
    #@-others
#@-node:gcross.20090827130017.1701:class TotalPotentialPathEnergyEstimates
#@-node:gcross.20090827130017.1692:Path estimates
#@-node:gcross.20090827130017.1678:Energy estimates
#@+node:gcross.20090827130017.1703:Position estimates
#@+node:gcross.20090827130017.1704:class AveragePositionEstimate
class AveragePositionEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    pass
#@-node:gcross.20090827130017.1704:class AveragePositionEstimate
#@+node:gcross.20090827130017.1705:class AverageAxialCoordinateEstimate
class AverageAxialCoordinateEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1706:__init__
    def __init__(self,axis,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        self.axis = axis
    #@-node:gcross.20090827130017.1706:__init__
    #@+node:gcross.20090827130017.1707:update
    def update(self):
        self.estimate += average(self.system.x[self.slice_number,:,self.axis])
    #@-node:gcross.20090827130017.1707:update
    #@-others
#@-node:gcross.20090827130017.1705:class AverageAxialCoordinateEstimate
#@+node:gcross.20090827130017.1708:class AverageAxialDistanceEstimate
class AverageAxialDistanceEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1709:__init__
    def __init__(self,axis,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        self.axis = axis
    #@-node:gcross.20090827130017.1709:__init__
    #@+node:gcross.20090827130017.1710:update
    def update(self):
        self.estimate += average(abs(self.system.x[self.slice_number,:,self.axis]))
    #@-node:gcross.20090827130017.1710:update
    #@-others
#@-node:gcross.20090827130017.1708:class AverageAxialDistanceEstimate
#@+node:gcross.20090827130017.1711:class AverageRadiusEstimate
class AverageRadiusEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1712:update
    def update(self):
        self.estimate += vpi.observables.compute_radius_average(self.system.x[self.slice_number])
    #@-node:gcross.20090827130017.1712:update
    #@-others
#@-node:gcross.20090827130017.1711:class AverageRadiusEstimate
#@+node:gcross.20090827130017.1713:class AveragePlaneRadiusEstimate
class AveragePlaneRadiusEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1714:__init__
    def __init__(self,plane_axis_1,plane_axis_2,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        assert plane_axis_1 >= 1
        assert plane_axis_2 >= 1
        assert not (plane_axis_1 == plane_axis_2)
        self.plane_axis_1 = plane_axis_1
        self.plane_axis_2 = plane_axis_2
    #@-node:gcross.20090827130017.1714:__init__
    #@+node:gcross.20090827130017.1715:update
    def update(self):
        self.estimate += vpi.observables.compute_plane_radius_average(self.system.x[self.slice_number],self.plane_axis_1,self.plane_axis_2)
    #@-node:gcross.20090827130017.1715:update
    #@-others
#@-node:gcross.20090827130017.1713:class AveragePlaneRadiusEstimate
#@+node:gcross.20090827130017.1716:class AverageRecipricalPlaneRadiusSquaredEstimate
class AverageRecipricalPlaneRadiusSquaredEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1717:__init__
    def __init__(self,plane_axis_1,plane_axis_2,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        assert plane_axis_1 >= 1
        assert plane_axis_2 >= 1
        assert not (plane_axis_1 == plane_axis_2)
        self.plane_axis_1 = plane_axis_1
        self.plane_axis_2 = plane_axis_2
    #@-node:gcross.20090827130017.1717:__init__
    #@+node:gcross.20090827130017.1718:update
    def update(self):
        self.estimate += vpi.observables.compute_recip_plane_r_sq_average(self.system.x,self.plane_axis_1,self.plane_axis_2)
    #@-node:gcross.20090827130017.1718:update
    #@-others
#@-node:gcross.20090827130017.1716:class AverageRecipricalPlaneRadiusSquaredEstimate
#@-node:gcross.20090827130017.1703:Position estimates
#@+node:gcross.20090827130017.1761:Rotation related estimates
#@+node:gcross.20090827130017.1719:class AverageAngularVelocityEstimate
class AverageAngularVelocityEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1720:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        self.estimate += average(first_derivatives)
    #@-node:gcross.20090827130017.1720:update
    #@-others
#@-node:gcross.20090827130017.1719:class AverageAngularVelocityEstimate
#@+node:gcross.20090827130017.1721:class AverageAngularVelocitySquaredEstimate
class AverageAngularVelocitySquaredEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1722:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        self.estimate += average(first_derivatives)
    #@-node:gcross.20090827130017.1722:update
    #@-others
#@-node:gcross.20090827130017.1721:class AverageAngularVelocitySquaredEstimate
#@+node:gcross.20090827130017.1723:class AverageRotationQuadraticTermEstimate
class AverageAngularVelocitySquaredEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1724:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        self.estimate += average(first_derivatives)
    #@-node:gcross.20090827130017.1724:update
    #@-others
#@-node:gcross.20090827130017.1723:class AverageRotationQuadraticTermEstimate
#@+node:gcross.20090827130017.1725:class StanderdDeviationAngularVelocityEstimate
class StandardDeviationAngularVelocityEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1726:update
    def update(self):
        system = self.system
        mean = float(system.number_of_rotating_particles) / system.number_of_particles
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        self.estimate += average((first_derivatives-mean)**2)
    #@-node:gcross.20090827130017.1726:update
    #@+node:gcross.20090827130017.1727:compute_total
    def compute_total(self):
        return sqrt(SingleAverageValueAtSliceEstimateAppendedToFile.compute_total(self))
    #@-node:gcross.20090827130017.1727:compute_total
    #@-others
#@-node:gcross.20090827130017.1725:class StanderdDeviationAngularVelocityEstimate
#@+node:gcross.20090827130017.1728:class AverageAngularSeparationEstimate
class AverageAngularSeparationEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1729:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        self.estimate += vpi.observables.compute_average_angular_separation(angles)
    #@-node:gcross.20090827130017.1729:update
    #@-others
#@-node:gcross.20090827130017.1728:class AverageAngularSeparationEstimate
#@+node:gcross.20090827130017.1730:class AverageNeighborAngularSeparationEstimate
class AverageNeighborAngularSeparationEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1731:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        self.estimate += vpi.observables.compute_avg_neighbor_angular_sep(angles)
    #@-node:gcross.20090827130017.1731:update
    #@-others
#@-node:gcross.20090827130017.1730:class AverageNeighborAngularSeparationEstimate
#@+node:gcross.20090827130017.1732:class AverageRotationQuadraticTermEstimate
class AverageRotationQuadraticTermEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1733:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            x,
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        term = (first_derivatives ** 2) / (x[:,system.rotation_plane_axis_2-1]**2+x[:,system.rotation_plane_axis_1-1]**2)
        self.estimate += average(term)
    #@-node:gcross.20090827130017.1733:update
    #@-others
#@-node:gcross.20090827130017.1732:class AverageRotationQuadraticTermEstimate
#@-node:gcross.20090827130017.1761:Rotation related estimates
#@-others
#@-node:gcross.20090827130017.1617:Observable classes
#@+node:gcross.20090827130017.1734:Functions
#@+node:gcross.20090827130017.1735:ensure_path_to_file_exists
def ensure_path_to_file_exists(path):
    directory, _ = os.path.split(path)
    if not os.path.exists(directory):
        os.makedirs(directory)
#@-node:gcross.20090827130017.1735:ensure_path_to_file_exists
#@-node:gcross.20090827130017.1734:Functions
#@+node:gcross.20090827130017.1736:class System
class System:
    #@    @+others
    #@+node:gcross.20090827130017.1737:Physics Functions
    #@+node:gcross.20090827130017.1738:compute_potential
    def compute_potential(self,x,_):
        x_sq = x**2
        U = array(dot(x_sq,self.harmonic_oscillator_coefficients)/2.0,dtype=double,order='Fortran')
        number_of_slices = x.shape[0]
        vpi.angular_momentum.accumulate_effective_potential2(
            x,self.lambda_,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            self.frame_angular_velocity,self.number_of_rotating_particles,
            U
        )
        gradU2  = zeros((number_of_slices,),dtype=double,order='Fortran')
        return U, gradU2, False
    #@-node:gcross.20090827130017.1738:compute_potential
    #@+node:gcross.20090827130017.1739:compute_trial
    def compute_trial(self,x,_):
        return -sum(dot(x**2,self.harmonic_oscillator_coefficients))/2
    #@-node:gcross.20090827130017.1739:compute_trial
    #@+node:gcross.20090827130017.1740:compute_trial_derivatives
    def compute_trial_derivatives(self,x,xij2):
        gradient_of_log_trial_fn = -x*self.harmonic_oscillator_coefficients
        laplacian_of_log_trial_fn = -sum(self.harmonic_oscillator_coefficients)*x.shape[0]
        return gradient_of_log_trial_fn, laplacian_of_log_trial_fn
    #@-node:gcross.20090827130017.1740:compute_trial_derivatives
    #@+node:gcross.20090828171041.1863:compute_greens_function
    def compute_greens_function(self,x,xij2,U,gradU2,lam,dt,slice_start,slice_end,particle_number):
        return vpi.gfn.gfn2_sp(slice_start,slice_end,particle_number,U,dt)
    #@-node:gcross.20090828171041.1863:compute_greens_function
    #@-node:gcross.20090827130017.1737:Physics Functions
    #@+node:gcross.20090827130017.1741:Observable management
    #@+node:gcross.20090827130017.1742:add_observable
    def add_observable(self,observable):
        self.observables.append(observable)
        observable.system = self
    #@-node:gcross.20090827130017.1742:add_observable
    #@+node:gcross.20090827130017.1743:total_and_write_observables
    def total_and_write_observables(self):
        for observable in self.observables:
            observable.total_and_write()
    #@-node:gcross.20090827130017.1743:total_and_write_observables
    #@-node:gcross.20090827130017.1741:Observable management
    #@+node:gcross.20090827130017.1744:__init__
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

        self.number_of_observations = self.total_number_of_observations // number_of_processors + 1
        self.total_number_of_observations = self.number_of_observations * number_of_processors

        self.number_of_thermalizations_per_observation = number_of_particles * number_of_slices // self.dM

        assert (number_of_slices % 2 == 0 and number_of_slices % 4 == 2)
        self.center_slice_number = number_of_slices // 2

        self.observables = []
    #@-node:gcross.20090827130017.1744:__init__
    #@+node:gcross.20090827130017.1745:run
    def run(self):
        #@    << Stash properties into local variables >>
        #@+node:gcross.20090827130017.1746:<< Stash properties into local variables >>
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
        compute_greens_function = self.compute_greens_function
        observables = self.observables
        #@-node:gcross.20090827130017.1746:<< Stash properties into local variables >>
        #@nl
        #@    << Prethermalize the system >>
        #@+node:gcross.20090827130017.1747:<< Prethermalize the system >>
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
            compute_potential,compute_trial,compute_greens_function
        )
        #@nonl
        #@-node:gcross.20090827130017.1747:<< Prethermalize the system >>
        #@nl
        #@    << Main iteration >>
        #@+node:gcross.20090827130017.1748:<< Main iteration >>
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
                compute_potential,compute_trial,compute_greens_function,
            )
            for observable in observables:
                observable.update()
            if (number_completed % decile == 0) and (my_rank == 0):
                print "{0:.0%} complete;  local bridge move acceptance rate = {1:.0%}, local rigid move acceptance rate = {2:.0%}".format(
                    float(number_completed)/self.number_of_observations,
                    float(move_type_accepted_counts[0])/move_type_attempted_counts[0],
                    float(move_type_accepted_counts[1])/move_type_attempted_counts[1],
                )
        #@-node:gcross.20090827130017.1748:<< Main iteration >>
        #@nl
    #@-node:gcross.20090827130017.1745:run
    #@-others
#@-node:gcross.20090827130017.1736:class System
#@-others
#@-node:gcross.20090827130017.1614:@thin system.py
#@-leo
