#@+leo-ver=4-thin
#@+node:gcross.20090828201103.1786:@thin system.py
#@<< Imports >>
#@+node:gcross.20090828201103.1787:<< Imports >>
import vpi

from numpy import *
from numpy.random import rand

import os
import os.path

import sys

import itertools
from itertools import izip
#@-node:gcross.20090828201103.1787:<< Imports >>
#@nl

#@<< MPI Initialization >>
#@+node:gcross.20090828201103.1788:<< MPI Initialization >>
from mpi4py import MPI

comm = MPI.COMM_WORLD

number_of_processors = comm.Get_size()
my_rank = comm.Get_rank()
#@-node:gcross.20090828201103.1788:<< MPI Initialization >>
#@nl

#@+others
#@+node:gcross.20090828201103.1789:Observable classes
#@+others
#@+node:gcross.20090828201103.1790:Base classes
#@+node:gcross.20090828201103.1791:class Observable
class Observable(object):
    #@    @+others
    #@+node:gcross.20090828201103.1792:total_and_write
    def total_and_write(self):
        totals = self.compute_total()
        if my_rank == 0:
            self.write_out_totals(totals)
    #@-node:gcross.20090828201103.1792:total_and_write
    #@-others
#@-node:gcross.20090828201103.1791:class Observable
#@+node:gcross.20090828201103.1793:class AverageValuesEstimate
class AverageValuesEstimate(Observable):
    #@    @+others
    #@+node:gcross.20090828201103.1794:(fields)
    total_number_of_observations = 0
    #@-node:gcross.20090828201103.1794:(fields)
    #@+node:gcross.20090828201103.1795:__init__
    def __init__(self,*shape):
        self.estimate = zeros(shape,dtype=double)
        self.estimate_squared = zeros(shape,dtype=double)
    #@-node:gcross.20090828201103.1795:__init__
    #@+node:gcross.20090828201103.1796:total_and_write
    def write_out_totals(self,totals_and_errors):
        totals, totals_squared = totals_and_errors
        errors = sqrt((totals_squared-totals**2)/(self.total_number_of_observations * number_of_processors))
        self.write_out(totals,errors)
    #@-node:gcross.20090828201103.1796:total_and_write
    #@+node:gcross.20090828201103.1797:compute_total
    def compute_total(self):
        total_estimate_and_square = zeros(2*prod(self.estimate.shape),dtype='d',order='Fortran')
        comm.Reduce((array([self.estimate,self.estimate_squared]).ravel(),MPI.DOUBLE),(total_estimate_and_square,MPI.DOUBLE))
        total_estimate_and_square /= (self.total_number_of_observations * number_of_processors)
        return total_estimate_and_square.reshape((2,)+self.estimate.shape)
    #@-node:gcross.20090828201103.1797:compute_total
    #@+node:gcross.20090828201103.1798:add
    def add(self,estimate):
        self.total_number_of_observations += 1
        self.estimate += estimate
        self.estimate_squared += estimate**2
    #@-node:gcross.20090828201103.1798:add
    #@-others
#@-node:gcross.20090828201103.1793:class AverageValuesEstimate
#@+node:gcross.20090828201103.1799:class SingleAverageValueEstimate
class SingleAverageValueEstimate(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1800:__init__
    def __init__(self):
        AverageValuesEstimate.__init__(self,1)
    #@-node:gcross.20090828201103.1800:__init__
    #@-others
#@-node:gcross.20090828201103.1799:class SingleAverageValueEstimate
#@+node:gcross.20090828201103.1801:class EstimatesAppendedToFile
class EstimatesAppendedToFile(Observable):
    #@    @+others
    #@+node:gcross.20090828201103.1802:__init__
    def __init__(self,filename,label):
        self.filename = filename
        self.label = label
    #@-node:gcross.20090828201103.1802:__init__
    #@+node:gcross.20090828201103.1803:write_out
    def write_out(self,total_estimate,total_estimate_error):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"a") as f:
            print >> f, self.label,
            for value, error in izip(total_estimate,total_estimate_error):
                print >> f, value, error,
            print >> f
    #@-node:gcross.20090828201103.1803:write_out
    #@-others
#@-node:gcross.20090828201103.1801:class EstimatesAppendedToFile
#@+node:gcross.20090828201103.1804:class SingleAverageValueEstimateAppendedToFile
class SingleAverageValueEstimateAppendedToFile(SingleAverageValueEstimate,EstimatesAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1805:__init__
    def __init__(self,filename,label):
        SingleAverageValueEstimate.__init__(self)
        EstimatesAppendedToFile.__init__(self,filename,label)
    #@-node:gcross.20090828201103.1805:__init__
    #@-others
#@-node:gcross.20090828201103.1804:class SingleAverageValueEstimateAppendedToFile
#@+node:gcross.20090828201103.1806:class SingleAverageValueAtSliceEstimateAppendedToFile
class SingleAverageValueAtSliceEstimateAppendedToFile(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1807:__init__
    def __init__(self,slice_number,label,filename):
        SingleAverageValueEstimateAppendedToFile.__init__(self,label,filename)
        self.slice_number = slice_number
    #@-node:gcross.20090828201103.1807:__init__
    #@-others
#@-node:gcross.20090828201103.1806:class SingleAverageValueAtSliceEstimateAppendedToFile
#@-node:gcross.20090828201103.1790:Base classes
#@+node:gcross.20090828201103.1808:Histograms
#@+node:gcross.20090828201103.1809:class Histogram
class Histogram(Observable):
    #@    @+others
    #@+node:gcross.20090828201103.1810:compute_total
    def compute_total(self):
        total_histogram = zeros(self.histogram.shape,dtype='i',order='Fortran')
        comm.Reduce((self.histogram,MPI.INT),(total_histogram,MPI.INT))
        return total_histogram
    #@-node:gcross.20090828201103.1810:compute_total
    #@+node:gcross.20090828201103.1811:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram))
        bin_width = float(self.right-self.left)/self.number_of_bins
        current = float(self.left)+bin_width/2
        with open(self.filename,"w") as f:
            for count in histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090828201103.1811:write_out_totals
    #@-others
#@-node:gcross.20090828201103.1809:class Histogram
#@+node:gcross.20090828201103.1812:class PositionDensity1DHistogram
class PositionDensity1DHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1813:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filenames):
        assert len(left) == len(right)
        self.left = array(left,dtype=double,order='Fortran')
        self.right = array(right,dtype=double,order='Fortran')
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((len(left),number_of_bins),dtype='i',order='Fortran')
        self.filenames = filenames
    #@-node:gcross.20090828201103.1813:__init__
    #@+node:gcross.20090828201103.1814:update
    def update(self):
        vpi.histograms.accumulate_1d_densities(
            self.system.x[self.slice_number],
            self.left,self.right,
            self.histogram
        )
    #@-node:gcross.20090828201103.1814:update
    #@+node:gcross.20090828201103.1815:write_out_totals
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
    #@-node:gcross.20090828201103.1815:write_out_totals
    #@-others
#@-node:gcross.20090828201103.1812:class PositionDensity1DHistogram
#@+node:gcross.20090828201103.1816:class RadialDensityHistogram
class RadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1817:(fields)
    left = 0
    #@-node:gcross.20090828201103.1817:(fields)
    #@+node:gcross.20090828201103.1818:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1818:__init__
    #@+node:gcross.20090828201103.1819:update
    def update(self):
        vpi.histograms.accumulate_radial_densities(
            self.system.x[self.slice_number],
            self.maximum_radius,
            self.histogram
        )
    #@-node:gcross.20090828201103.1819:update
    #@-others
#@-node:gcross.20090828201103.1816:class RadialDensityHistogram
#@+node:gcross.20090828201103.1820:class PlaneRadialDensityHistogram
class PlaneRadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1821:(fields)
    left = 0
    #@-node:gcross.20090828201103.1821:(fields)
    #@+node:gcross.20090828201103.1822:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1822:__init__
    #@+node:gcross.20090828201103.1823:update
    def update(self):
        vpi.histograms.accumulate_plane_radial_densities(
            self.system.x[self.slice_number],
            self.maximum_radius,
            self.system.rotation_plane_axis_1,
            self.system.rotation_plane_axis_2,
            self.histogram
        )
    #@-node:gcross.20090828201103.1823:update
    #@-others
#@-node:gcross.20090828201103.1820:class PlaneRadialDensityHistogram
#@+node:gcross.20090828201103.1824:class RecipricalRadiusSquaredDensityHistogram
class RecipricalRadiusSquaredDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1825:(fields)
    left = 0
    #@-node:gcross.20090828201103.1825:(fields)
    #@+node:gcross.20090828201103.1826:__init__
    def __init__(self,slice_number,maximum_value,number_of_bins,filename):
        self.right = self.maximum_value = maximum_value
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1826:__init__
    #@+node:gcross.20090828201103.1827:update
    def update(self):
        vpi.histograms.accumulate_reciprical_radius_squared_densities(
            self.system.x[self.slice_number],
            self.maximum_value,
            self.histogram
        )
    #@-node:gcross.20090828201103.1827:update
    #@-others
#@-node:gcross.20090828201103.1824:class RecipricalRadiusSquaredDensityHistogram
#@+node:gcross.20090828201103.1828:class RecipricalPlaneRadiusSquaredDensityHistogram
class RecipricalPlaneRadiusSquaredDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1829:(fields)
    left = 0
    #@-node:gcross.20090828201103.1829:(fields)
    #@+node:gcross.20090828201103.1830:__init__
    def __init__(self,slice_number,maximum_value,number_of_bins,filename):
        self.right = self.maximum_value = maximum_value
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1830:__init__
    #@+node:gcross.20090828201103.1831:update
    def update(self):
        vpi.histograms.accumulate_recip_plane_r_sq_densities(
            self.system.x[self.slice_number],
            self.maximum_value,
            self.system.rotation_plane_axis_1,
            self.system.rotation_plane_axis_2,
            self.histogram
        )
    #@-node:gcross.20090828201103.1831:update
    #@+node:gcross.20090828201103.1832:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(self.maximum_value)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090828201103.1832:write_out_totals
    #@-others
#@-node:gcross.20090828201103.1828:class RecipricalPlaneRadiusSquaredDensityHistogram
#@+node:gcross.20090828201103.1833:class AngularSeparationDensityHistogram
class AngularSeparationDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1834:(fields)
    left = 0
    right = 2*pi
    #@-node:gcross.20090828201103.1834:(fields)
    #@+node:gcross.20090828201103.1835:__init__
    def __init__(self,slice_number,number_of_bins,filename):
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1835:__init__
    #@+node:gcross.20090828201103.1836:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        vpi.histograms.accumulate_angular_separation_densities(
            angles,
            self.histogram
        )
    #@-node:gcross.20090828201103.1836:update
    #@+node:gcross.20090828201103.1837:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(2*pi)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090828201103.1837:write_out_totals
    #@-others
#@-node:gcross.20090828201103.1833:class AngularSeparationDensityHistogram
#@+node:gcross.20090828201103.1838:class NeighborAngularSeparationDensityHistogram
class NeighborAngularSeparationDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1839:(fields)
    left = 0
    right = 2*pi
    #@-node:gcross.20090828201103.1839:(fields)
    #@+node:gcross.20090828201103.1840:__init__
    def __init__(self,slice_number,number_of_bins,filename):
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1840:__init__
    #@+node:gcross.20090828201103.1841:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        vpi.histograms.accumulate_neighbor_angular_separation_densities(
            angles,
            self.histogram
        )
    #@-node:gcross.20090828201103.1841:update
    #@-others
#@-node:gcross.20090828201103.1838:class NeighborAngularSeparationDensityHistogram
#@+node:gcross.20090828201103.1842:class AngularVelocityHistogram
class AngularVelocityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1843:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1843:__init__
    #@+node:gcross.20090828201103.1844:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        vpi.histograms.accumulate(first_derivatives,self.left,self.right,self.histogram)
    #@-node:gcross.20090828201103.1844:update
    #@-others
#@-node:gcross.20090828201103.1842:class AngularVelocityHistogram
#@+node:gcross.20090828201103.1845:class AngularVelocitySquaredHistogram
class AngularVelocitySquaredHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1846:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1846:__init__
    #@+node:gcross.20090828201103.1847:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        vpi.histograms.accumulate(first_derivatives,self.left,self.right,self.histogram)
    #@-node:gcross.20090828201103.1847:update
    #@-others
#@-node:gcross.20090828201103.1845:class AngularVelocitySquaredHistogram
#@+node:gcross.20090828201103.1848:class RotationQuadraticTermHistogram
class RotationQuadraticTermHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1849:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1849:__init__
    #@+node:gcross.20090828201103.1850:update
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
    #@-node:gcross.20090828201103.1850:update
    #@-others
#@-node:gcross.20090828201103.1848:class RotationQuadraticTermHistogram
#@+node:gcross.20090828201103.1851:class AngularVelocityAndRadiusHistogram
class AngularVelocityAndRadiusHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090828201103.1852:__init__
    def __init__(self,slice_number,maximum_angular_velocity,maximum_radius,number_of_bins,filename):
        self.maximum_angular_velocity = maximum_angular_velocity
        self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,)*2,dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090828201103.1852:__init__
    #@+node:gcross.20090828201103.1853:update
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
    #@-node:gcross.20090828201103.1853:update
    #@+node:gcross.20090828201103.1854:write_out_totals
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
    #@-node:gcross.20090828201103.1854:write_out_totals
    #@-others
#@-node:gcross.20090828201103.1851:class AngularVelocityAndRadiusHistogram
#@+node:gcross.20090830224709.2057:class ParticleSeparationHistogram
class ParticleSeparationHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090830224709.2058:(fields)
    left = 0
    #@-node:gcross.20090830224709.2058:(fields)
    #@+node:gcross.20090830224709.2059:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_value = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090830224709.2059:__init__
    #@+node:gcross.20090830224709.2060:update
    def update(self):
        vpi.histograms.accumulate_particle_separation_densities(
            self.system.xij2[self.slice_number],
            self.maximum_value,
            self.histogram
        )
    #@-node:gcross.20090830224709.2060:update
    #@-others
#@-node:gcross.20090830224709.2057:class ParticleSeparationHistogram
#@-node:gcross.20090828201103.1808:Histograms
#@+node:gcross.20090828201103.1855:Energy estimates
#@+node:gcross.20090828201103.1856:class TotalEnergyEstimate
class TotalEnergyEstimate(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1857:update
    def update(self):
        system = self.system
        for slice_number in [0,-1]:
            gradient_of_log_trial_fn, laplacian_of_log_trial_fn = system.compute_trial_derivatives(system.x[slice_number],system.xij2[slice_number])
            self.add(
                vpi.observables.compute_local_energy_estimate(
                    system.U[slice_number],
                    gradient_of_log_trial_fn, laplacian_of_log_trial_fn,
                    system.lambda_,
                )
            )
    #@-node:gcross.20090828201103.1857:update
    #@-others
#@-node:gcross.20090828201103.1856:class TotalEnergyEstimate
#@+node:gcross.20090828201103.1858:Slice estimates
#@+node:gcross.20090828201103.1859:class SliceEnergyEstimate
class SliceEnergyEstimate(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1860:__init__
    def __init__(self,slice_number,filename,label):
        SingleAverageValueEstimateAppendedToFile.__init__(self,filename,label)
        self.slice_number = slice_number
    #@-node:gcross.20090828201103.1860:__init__
    #@-others
#@-node:gcross.20090828201103.1859:class SliceEnergyEstimate
#@+node:gcross.20090828201103.1861:class EffectivePotentialSliceEnergyEstimate
class EffectivePotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1862:update
    def update(self):
        self.add(self.system.compute_effective_potential(self.slice_number))
    #@-node:gcross.20090828201103.1862:update
    #@-others
#@-node:gcross.20090828201103.1861:class EffectivePotentialSliceEnergyEstimate
#@+node:gcross.20090828201103.1863:class PhysicalPotentialSliceEnergyEstimate
class PhysicalPotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1864:update
    def update(self):
            self.add(self.system.compute_physical_potential(self.slice_number))
    #@-node:gcross.20090828201103.1864:update
    #@-others
#@-node:gcross.20090828201103.1863:class PhysicalPotentialSliceEnergyEstimate
#@+node:gcross.20090828201103.1865:class TotalPotentialSliceEnergyEstimate
class TotalPotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1866:update
    def update(self):
        self.add(self.system.compute_total_potential(self.slice_number))
    #@-node:gcross.20090828201103.1866:update
    #@-others
#@-node:gcross.20090828201103.1865:class TotalPotentialSliceEnergyEstimate
#@-node:gcross.20090828201103.1858:Slice estimates
#@+node:gcross.20090828201103.1867:Path estimates
#@+node:gcross.20090828201103.1868:class PathEnergyEstimates
class PathEnergyEstimates(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1869:(fields)
    estimates = 0
    #@-node:gcross.20090828201103.1869:(fields)
    #@+node:gcross.20090828201103.1870:__init__
    def __init__(self,filename):
        self.filename = filename
    #@-node:gcross.20090828201103.1870:__init__
    #@+node:gcross.20090828201103.1871:write_out_totals
    def write_out_totals(self,total_estimates):
        ensure_path_to_file_exists(self.filename)
        center_slice_number = self.system.center_slice_number
        with open(self.filename,"w") as f:
            for slice_number, estimate in enumerate(total_estimates):
                print >> f, center_slice_number-abs(center_slice_number-slice_number), estimate, slice_number
    #@-node:gcross.20090828201103.1871:write_out_totals
    #@-others
#@-node:gcross.20090828201103.1868:class PathEnergyEstimates
#@+node:gcross.20090828201103.1872:class EffectivePotentialPathEnergyEstimates
class EffectivePotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090828201103.1873:update
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
    #@-node:gcross.20090828201103.1873:update
    #@-others
#@-node:gcross.20090828201103.1872:class EffectivePotentialPathEnergyEstimates
#@+node:gcross.20090828201103.1874:class PhysicalPotentialPathEnergyEstimates
class PhysicalPotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090828201103.1875:update
    def update(self):
        self.estimates += sum(dot(self.system.x**2,self.system.harmonic_oscillator_coefficients),axis=-1)/2.0
    #@-node:gcross.20090828201103.1875:update
    #@-others
#@-node:gcross.20090828201103.1874:class PhysicalPotentialPathEnergyEstimates
#@+node:gcross.20090828201103.1876:class TotalPotentialPathEnergyEstimates
class TotalPotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090828201103.1877:update
    def update(self):
        self.estimates += sum(self.system.U,axis=-1)
    #@-node:gcross.20090828201103.1877:update
    #@-others
#@-node:gcross.20090828201103.1876:class TotalPotentialPathEnergyEstimates
#@-node:gcross.20090828201103.1867:Path estimates
#@-node:gcross.20090828201103.1855:Energy estimates
#@+node:gcross.20090828201103.1878:Position estimates
#@+node:gcross.20090828201103.1879:class AveragePositionEstimate
class AveragePositionEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    pass
#@-node:gcross.20090828201103.1879:class AveragePositionEstimate
#@+node:gcross.20090828201103.1880:class AverageAxialCoordinateEstimate
class AverageAxialCoordinateEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1881:__init__
    def __init__(self,axis,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        self.axis = axis
    #@-node:gcross.20090828201103.1881:__init__
    #@+node:gcross.20090828201103.1882:update
    def update(self):
        self.add(average(self.system.x[self.slice_number,:,self.axis]))
    #@-node:gcross.20090828201103.1882:update
    #@-others
#@-node:gcross.20090828201103.1880:class AverageAxialCoordinateEstimate
#@+node:gcross.20090828201103.1883:class AverageAxialDistanceEstimate
class AverageAxialDistanceEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1884:__init__
    def __init__(self,axis,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        self.axis = axis
    #@-node:gcross.20090828201103.1884:__init__
    #@+node:gcross.20090828201103.1885:update
    def update(self):
        self.add(average(abs(self.system.x[self.slice_number,:,self.axis])))
    #@-node:gcross.20090828201103.1885:update
    #@-others
#@-node:gcross.20090828201103.1883:class AverageAxialDistanceEstimate
#@+node:gcross.20090828201103.1886:class AverageRadiusEstimate
class AverageRadiusEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1887:update
    def update(self):
        self.add(vpi.observables.compute_radius_average(self.system.x[self.slice_number]))
    #@-node:gcross.20090828201103.1887:update
    #@-others
#@-node:gcross.20090828201103.1886:class AverageRadiusEstimate
#@+node:gcross.20090828201103.1888:class AveragePlaneRadiusEstimate
class AveragePlaneRadiusEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1889:__init__
    def __init__(self,plane_axis_1,plane_axis_2,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        assert plane_axis_1 >= 1
        assert plane_axis_2 >= 1
        assert not (plane_axis_1 == plane_axis_2)
        self.plane_axis_1 = plane_axis_1
        self.plane_axis_2 = plane_axis_2
    #@-node:gcross.20090828201103.1889:__init__
    #@+node:gcross.20090828201103.1890:update
    def update(self):
        self.add( vpi.observables.compute_plane_radius_average(self.system.x[self.slice_number],self.plane_axis_1,self.plane_axis_2))
    #@-node:gcross.20090828201103.1890:update
    #@-others
#@-node:gcross.20090828201103.1888:class AveragePlaneRadiusEstimate
#@+node:gcross.20090828201103.1891:class AverageRecipricalPlaneRadiusSquaredEstimate
class AverageRecipricalPlaneRadiusSquaredEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090828201103.1892:__init__
    def __init__(self,plane_axis_1,plane_axis_2,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        assert plane_axis_1 >= 1
        assert plane_axis_2 >= 1
        assert not (plane_axis_1 == plane_axis_2)
        self.plane_axis_1 = plane_axis_1
        self.plane_axis_2 = plane_axis_2
    #@-node:gcross.20090828201103.1892:__init__
    #@+node:gcross.20090828201103.1893:update
    def update(self):
        self.add(vpi.observables.compute_recip_plane_r_sq_average(self.system.x,self.plane_axis_1,self.plane_axis_2))
    #@-node:gcross.20090828201103.1893:update
    #@-others
#@-node:gcross.20090828201103.1891:class AverageRecipricalPlaneRadiusSquaredEstimate
#@+node:gcross.20090830224709.2064:class AverageParticleSeparationEstimate
class AverageParticleSeparationEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090830224709.2065:__init__
    def __init__(self,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
    #@-node:gcross.20090830224709.2065:__init__
    #@+node:gcross.20090830224709.2066:update
    def update(self):
        self.add(vpi.observables.compute_particle_separation_average(self.system.xij2[self.slice_number]))
    #@-node:gcross.20090830224709.2066:update
    #@-others
#@-node:gcross.20090830224709.2064:class AverageParticleSeparationEstimate
#@-node:gcross.20090828201103.1878:Position estimates
#@+node:gcross.20090828201103.1894:Rotation related estimates
#@+node:gcross.20090828201103.1895:class AverageAngularVelocityEstimate
class AverageAngularVelocityEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1896:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        self.add(average(first_derivatives))
    #@-node:gcross.20090828201103.1896:update
    #@-others
#@-node:gcross.20090828201103.1895:class AverageAngularVelocityEstimate
#@+node:gcross.20090828201103.1897:class AverageAngularVelocitySquaredEstimate
class AverageAngularVelocitySquaredEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1898:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        self.add(average(first_derivatives))
    #@-node:gcross.20090828201103.1898:update
    #@-others
#@-node:gcross.20090828201103.1897:class AverageAngularVelocitySquaredEstimate
#@+node:gcross.20090828201103.1899:class AverageRotationQuadraticTermEstimate
class AverageAngularVelocitySquaredEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1900:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        self.add(average(first_derivatives))
    #@-node:gcross.20090828201103.1900:update
    #@-others
#@-node:gcross.20090828201103.1899:class AverageRotationQuadraticTermEstimate
#@+node:gcross.20090828201103.1901:class AverageAngularSeparationEstimate
class AverageAngularSeparationEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1902:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        self.add(vpi.observables.compute_average_angular_separation(angles))
    #@-node:gcross.20090828201103.1902:update
    #@-others
#@-node:gcross.20090828201103.1901:class AverageAngularSeparationEstimate
#@+node:gcross.20090828201103.1903:class AverageNeighborAngularSeparationEstimate
class AverageNeighborAngularSeparationEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1904:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        self.add(vpi.observables.compute_avg_neighbor_angular_sep(angles))
    #@-node:gcross.20090828201103.1904:update
    #@-others
#@-node:gcross.20090828201103.1903:class AverageNeighborAngularSeparationEstimate
#@+node:gcross.20090828201103.1905:class AverageRotationQuadraticTermEstimate
class AverageRotationQuadraticTermEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090828201103.1906:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            x,
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        term = (first_derivatives ** 2) / (x[:,system.rotation_plane_axis_2-1]**2+x[:,system.rotation_plane_axis_1-1]**2)
        self.add(average(term))
    #@-node:gcross.20090828201103.1906:update
    #@-others
#@-node:gcross.20090828201103.1905:class AverageRotationQuadraticTermEstimate
#@-node:gcross.20090828201103.1894:Rotation related estimates
#@-others
#@-node:gcross.20090828201103.1789:Observable classes
#@+node:gcross.20090828201103.1907:Functions
#@+node:gcross.20090828201103.1908:ensure_path_to_file_exists
def ensure_path_to_file_exists(path):
    directory, _ = os.path.split(path)
    if not os.path.exists(directory):
        os.makedirs(directory)
#@-node:gcross.20090828201103.1908:ensure_path_to_file_exists
#@-node:gcross.20090828201103.1907:Functions
#@+node:gcross.20090828201103.1909:class System
class System:
    #@    @+others
    #@+node:gcross.20090828201103.1910:Physics Functions
    #@+node:gcross.20090901084550.2198:Potential
    #@+node:gcross.20090828201103.1911:compute_potential
    def compute_potential(self,x,xij2):
        U = zeros(x.shape[:2],dtype=double,order='Fortran')
        gradU2 = zeros(x.shape[:1],dtype=double,order='Fortran')
        if (vpi.xij.hard_sphere_violation(xij2,self.hard_sphere_radius_squared)):
            return U, gradU2, True
        else:
            vpi.harmonic_oscillator_3d.accumulate_potential(x,self.harmonic_oscillator_coefficients,U)
            self.accumulate_effective_potential(x,U)
            return U, gradU2, False
    #@-node:gcross.20090828201103.1911:compute_potential
    #@+node:gcross.20090901084550.2194:accumulate_effective_potential
    def accumulate_effective_potential(self,x,U):
        vpi.angular_momentum.accumulate_effective_potential(
            x,self.lambda_,
            self.rotation_plane_axis_1,self.rotation_plane_axis_2,
            self.frame_angular_velocity,self.number_of_rotating_particles,
            U
        )
    #@-node:gcross.20090901084550.2194:accumulate_effective_potential
    #@+node:gcross.20090901084550.2195:compute_effective_potential
    def compute_effective_potential(self,slice_number):
        return self.compute_total_potential(slice_number)-self.compute_physical_potential(slice_number)
    #@-node:gcross.20090901084550.2195:compute_effective_potential
    #@+node:gcross.20090901084550.2197:compute_physical_potential
    def compute_physical_potential(self,slice_number):
        x_sq = self.x[slice_number]**2
        U = array(dot(x_sq,self.harmonic_oscillator_coefficients)/2.0,dtype=double,order='Fortran')
        return sum(U)
    #@-node:gcross.20090901084550.2197:compute_physical_potential
    #@+node:gcross.20090901084550.2200:compute_total_potential
    def compute_total_potential(self,slice_number):
        return sum(self.U[slice_number])
    #@-node:gcross.20090901084550.2200:compute_total_potential
    #@-node:gcross.20090901084550.2198:Potential
    #@+node:gcross.20090901084550.2201:Trial
    #@+node:gcross.20090828201103.1912:compute_trial
    def compute_trial(self,x,xij2):
        weight1, = vpi.harmonic_oscillator_3d.compute_trial_weight(x,self.harmonic_oscillator_coefficients),
        weight2, reject2 = vpi.hard_sphere_interaction.compute_trial_weight(xij2,self.hard_sphere_radius)
        return (weight1+weight2), reject2
    #@-node:gcross.20090828201103.1912:compute_trial
    #@+node:gcross.20090828201103.1913:compute_trial_derivatives
    def compute_trial_derivatives(self,x,xij2):
        gradient_of_log_trial_fn = zeros(x.shape,dtype=double,order='Fortran')
        laplacian_of_log_trial_fn = zeros((),dtype=double,order='Fortran')
        vpi.harmonic_oscillator_3d.accumulate_trial_derivatives(
            x,self.harmonic_oscillator_coefficients,
            gradient_of_log_trial_fn,laplacian_of_log_trial_fn
        )
        vpi.hard_sphere_interaction.accumulate_trial_derivatives(
            x,xij2,self.hard_sphere_radius,
            gradient_of_log_trial_fn,laplacian_of_log_trial_fn
        )
        return gradient_of_log_trial_fn, laplacian_of_log_trial_fn
    #@-node:gcross.20090828201103.1913:compute_trial_derivatives
    #@-node:gcross.20090901084550.2201:Trial
    #@+node:gcross.20090828201103.1914:compute_greens_function
    def compute_greens_function(self,x,xij2,U,gradU2,lam,dt,slice_start,slice_end,particle_number):
        return (
            vpi.gfn.gfn2_sp(slice_start,slice_end,particle_number,U,dt)
         +  vpi.gfn.gfn_hard_sphere_contribution(xij2,dt,self.hard_sphere_radius,slice_start,slice_end,particle_number)
         )
    #@-node:gcross.20090828201103.1914:compute_greens_function
    #@-node:gcross.20090828201103.1910:Physics Functions
    #@+node:gcross.20090828201103.1915:Observable management
    #@+node:gcross.20090828201103.1916:add_observable
    def add_observable(self,observable):
        self.observables.append(observable)
        observable.system = self
    #@-node:gcross.20090828201103.1916:add_observable
    #@+node:gcross.20090828201103.1917:total_and_write_observables
    def total_and_write_observables(self):
        for observable in self.observables:
            observable.total_and_write()
    #@-node:gcross.20090828201103.1917:total_and_write_observables
    #@-node:gcross.20090828201103.1915:Observable management
    #@+node:gcross.20090828201103.1918:__init__
    def __init__(self,**keywords):
        self.__dict__.update(keywords)

        number_of_slices = self.number_of_slices
        number_of_particles = self.number_of_particles
        number_of_dimensions = self.number_of_dimensions

        vpi.rand_utils.init_seed(my_rank)

        hard_sphere_condition_violated = True
        self.hard_sphere_radius_squared = self.hard_sphere_radius**2
        number_of_attempts = 0
        while(hard_sphere_condition_violated and number_of_attempts < 100):
            number_of_attempts += 1
            self.x = vpi.lattice.make_lattice(self.initial_particle_distribution_size,number_of_slices,number_of_particles,number_of_dimensions)
            self.xij2 = zeros((number_of_slices,number_of_particles,number_of_particles),dtype=double,order='Fortran')
            vpi.xij.update_xij(self.xij2,self.x)
            hard_sphere_condition_violated = vpi.xij.hard_sphere_violation(self.xij2[0:1],self.hard_sphere_radius_squared)
        if(number_of_attempts == 100):
            print >> sys.stderr, "Failed in 100 attempts to construct a system where no two particles violated the hard sphere condition."
            comm.Abort(-1)

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
    #@-node:gcross.20090828201103.1918:__init__
    #@+node:gcross.20090828201103.1919:run
    def run(self):
        #@    << Stash properties into local variables >>
        #@+node:gcross.20090828201103.1920:<< Stash properties into local variables >>
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
        #@-node:gcross.20090828201103.1920:<< Stash properties into local variables >>
        #@nl
        #@    << Prethermalize the system >>
        #@+node:gcross.20090828201103.1921:<< Prethermalize the system >>
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
        #@-node:gcross.20090828201103.1921:<< Prethermalize the system >>
        #@nl
        #@    << Main iteration >>
        #@+node:gcross.20090828201103.1922:<< Main iteration >>
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
        #@-node:gcross.20090828201103.1922:<< Main iteration >>
        #@nl
    #@-node:gcross.20090828201103.1919:run
    #@-others
#@-node:gcross.20090828201103.1909:class System
#@-others
#@-node:gcross.20090828201103.1786:@thin system.py
#@-leo
