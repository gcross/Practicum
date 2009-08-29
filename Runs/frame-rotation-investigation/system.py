#@+leo-ver=4-thin
#@+node:gcross.20090827130017.1614:@thin system.py
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
#@+node:gcross.20090827130017.1895:Observable classes
#@+others
#@+node:gcross.20090827130017.1896:Base classes
#@+node:gcross.20090827130017.1897:class Observable
class Observable(object):
    #@    @+others
    #@+node:gcross.20090827130017.2016:total_and_write
    def total_and_write(self):
        totals = self.compute_total()
        if my_rank == 0:
            self.write_out_totals(totals)
    #@-node:gcross.20090827130017.2016:total_and_write
    #@-others
#@-node:gcross.20090827130017.1897:class Observable
#@+node:gcross.20090827130017.1899:class AverageValuesEstimate
class AverageValuesEstimate(Observable):
    #@    @+others
    #@+node:gcross.20090828095451.2961:(fields)
    total_number_of_observations = 0
    #@-node:gcross.20090828095451.2961:(fields)
    #@+node:gcross.20090827130017.2017:__init__
    def __init__(self,*shape):
        self.estimate = zeros(shape,dtype=double)
        self.estimate_squared = zeros(shape,dtype=double)
    #@-node:gcross.20090827130017.2017:__init__
    #@+node:gcross.20090827130017.1898:total_and_write
    def write_out_totals(self,totals_and_errors):
        totals, totals_squared = totals_and_errors
        errors = sqrt((totals_squared-totals**2)/(self.total_number_of_observations * number_of_processors))
        self.write_out(totals,errors)
    #@-node:gcross.20090827130017.1898:total_and_write
    #@+node:gcross.20090827130017.1900:compute_total
    def compute_total(self):
        total_estimate_and_square = zeros(2*prod(self.estimate.shape),dtype='d',order='Fortran')
        comm.Reduce((array([self.estimate,self.estimate_squared]).ravel(),MPI.DOUBLE),(total_estimate_and_square,MPI.DOUBLE))
        total_estimate_and_square /= (self.total_number_of_observations * number_of_processors)
        return total_estimate_and_square.reshape((2,)+self.estimate.shape)
    #@-node:gcross.20090827130017.1900:compute_total
    #@+node:gcross.20090827130017.2014:add
    def add(self,estimate):
        self.total_number_of_observations += 1
        self.estimate += estimate
        self.estimate_squared += estimate**2
    #@-node:gcross.20090827130017.2014:add
    #@-others
#@-node:gcross.20090827130017.1899:class AverageValuesEstimate
#@+node:gcross.20090827130017.1901:class SingleAverageValueEstimate
class SingleAverageValueEstimate(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1902:__init__
    def __init__(self):
        AverageValuesEstimate.__init__(self,1)
    #@-node:gcross.20090827130017.1902:__init__
    #@-others
#@-node:gcross.20090827130017.1901:class SingleAverageValueEstimate
#@+node:gcross.20090827130017.1903:class EstimatesAppendedToFile
class EstimatesAppendedToFile(Observable):
    #@    @+others
    #@+node:gcross.20090827130017.1904:__init__
    def __init__(self,filename,label):
        self.filename = filename
        self.label = label
    #@-node:gcross.20090827130017.1904:__init__
    #@+node:gcross.20090827130017.1905:write_out
    def write_out(self,total_estimate,total_estimate_error):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"a") as f:
            print >> f, self.label,
            for value, error in izip(total_estimate,total_estimate_error):
                print >> f, value, error
            print >> f
    #@-node:gcross.20090827130017.1905:write_out
    #@-others
#@-node:gcross.20090827130017.1903:class EstimatesAppendedToFile
#@+node:gcross.20090827130017.1906:class SingleAverageValueEstimateAppendedToFile
class SingleAverageValueEstimateAppendedToFile(SingleAverageValueEstimate,EstimatesAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1907:__init__
    def __init__(self,filename,label):
        SingleAverageValueEstimate.__init__(self)
        EstimatesAppendedToFile.__init__(self,filename,label)
    #@-node:gcross.20090827130017.1907:__init__
    #@-others
#@-node:gcross.20090827130017.1906:class SingleAverageValueEstimateAppendedToFile
#@+node:gcross.20090827130017.1908:class SingleAverageValueAtSliceEstimateAppendedToFile
class SingleAverageValueAtSliceEstimateAppendedToFile(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1909:__init__
    def __init__(self,slice_number,label,filename):
        SingleAverageValueEstimateAppendedToFile.__init__(self,label,filename)
        self.slice_number = slice_number
    #@-node:gcross.20090827130017.1909:__init__
    #@-others
#@-node:gcross.20090827130017.1908:class SingleAverageValueAtSliceEstimateAppendedToFile
#@-node:gcross.20090827130017.1896:Base classes
#@+node:gcross.20090827130017.1910:Histograms
#@+node:gcross.20090827130017.1911:class Histogram
class Histogram(Observable):
    #@    @+others
    #@+node:gcross.20090827130017.1912:compute_total
    def compute_total(self):
        total_histogram = zeros(self.histogram.shape,dtype='i',order='Fortran')
        comm.Reduce((self.histogram,MPI.INT),(total_histogram,MPI.INT))
        return total_histogram
    #@-node:gcross.20090827130017.1912:compute_total
    #@+node:gcross.20090827130017.1913:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram))
        bin_width = float(self.right-self.left)/self.number_of_bins
        current = float(self.left)+bin_width/2
        with open(self.filename,"w") as f:
            for count in histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090827130017.1913:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1911:class Histogram
#@+node:gcross.20090827130017.1914:class PositionDensity1DHistogram
class PositionDensity1DHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1915:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filenames):
        assert len(left) == len(right)
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((len(left),number_of_bins),dtype='i',order='Fortran')
        self.filenames = filenames
    #@-node:gcross.20090827130017.1915:__init__
    #@+node:gcross.20090827130017.1916:update
    def update(self):
        vpi.histograms.accumulate_1d_densities(
            self.system.x[self.slice_number],
            self.left,self.right,
            self.histogram
        )
    #@-node:gcross.20090827130017.1916:update
    #@+node:gcross.20090827130017.1917:write_out_totals
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
    #@-node:gcross.20090827130017.1917:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1914:class PositionDensity1DHistogram
#@+node:gcross.20090827130017.1918:class RadialDensityHistogram
class RadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1919:(fields)
    left = 0
    #@-node:gcross.20090827130017.1919:(fields)
    #@+node:gcross.20090827130017.1920:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1920:__init__
    #@+node:gcross.20090827130017.1921:update
    def update(self):
        vpi.histograms.accumulate_radial_densities(
            self.system.x[self.slice_number],
            self.maximum_radius,
            self.histogram
        )
    #@-node:gcross.20090827130017.1921:update
    #@-others
#@-node:gcross.20090827130017.1918:class RadialDensityHistogram
#@+node:gcross.20090827130017.1922:class PlaneRadialDensityHistogram
class PlaneRadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1923:(fields)
    left = 0
    #@-node:gcross.20090827130017.1923:(fields)
    #@+node:gcross.20090827130017.1924:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1924:__init__
    #@+node:gcross.20090827130017.1925:update
    def update(self):
        vpi.histograms.accumulate_plane_radial_densities(
            self.system.x[self.slice_number],
            self.maximum_radius,
            self.system.rotation_plane_axis_1,
            self.system.rotation_plane_axis_2,
            self.histogram
        )
    #@-node:gcross.20090827130017.1925:update
    #@-others
#@-node:gcross.20090827130017.1922:class PlaneRadialDensityHistogram
#@+node:gcross.20090827130017.1926:class RecipricalRadiusSquaredDensityHistogram
class RecipricalRadiusSquaredDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1927:(fields)
    left = 0
    #@-node:gcross.20090827130017.1927:(fields)
    #@+node:gcross.20090827130017.1928:__init__
    def __init__(self,slice_number,maximum_value,number_of_bins,filename):
        self.right = self.maximum_value = maximum_value
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1928:__init__
    #@+node:gcross.20090827130017.1929:update
    def update(self):
        vpi.histograms.accumulate_reciprical_radius_squared_densities(
            self.system.x[self.slice_number],
            self.maximum_value,
            self.histogram
        )
    #@-node:gcross.20090827130017.1929:update
    #@-others
#@-node:gcross.20090827130017.1926:class RecipricalRadiusSquaredDensityHistogram
#@+node:gcross.20090827130017.1930:class RecipricalPlaneRadiusSquaredDensityHistogram
class RecipricalPlaneRadiusSquaredDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1931:(fields)
    left = 0
    #@-node:gcross.20090827130017.1931:(fields)
    #@+node:gcross.20090827130017.1932:__init__
    def __init__(self,slice_number,maximum_value,number_of_bins,filename):
        self.right = self.maximum_value = maximum_value
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1932:__init__
    #@+node:gcross.20090827130017.1933:update
    def update(self):
        vpi.histograms.accumulate_recip_plane_r_sq_densities(
            self.system.x[self.slice_number],
            self.maximum_value,
            self.system.rotation_plane_axis_1,
            self.system.rotation_plane_axis_2,
            self.histogram
        )
    #@-node:gcross.20090827130017.1933:update
    #@+node:gcross.20090827130017.1934:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(self.maximum_value)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090827130017.1934:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1930:class RecipricalPlaneRadiusSquaredDensityHistogram
#@+node:gcross.20090827130017.1935:class AngularSeparationDensityHistogram
class AngularSeparationDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1936:(fields)
    left = 0
    right = 2*pi
    #@-node:gcross.20090827130017.1936:(fields)
    #@+node:gcross.20090827130017.1937:__init__
    def __init__(self,slice_number,number_of_bins,filename):
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1937:__init__
    #@+node:gcross.20090827130017.1938:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        vpi.histograms.accumulate_angular_separation_densities(
            angles,
            self.histogram
        )
    #@-node:gcross.20090827130017.1938:update
    #@+node:gcross.20090827130017.1939:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(2*pi)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090827130017.1939:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1935:class AngularSeparationDensityHistogram
#@+node:gcross.20090827130017.1940:class NeighborAngularSeparationDensityHistogram
class NeighborAngularSeparationDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1941:(fields)
    left = 0
    right = 2*pi
    #@-node:gcross.20090827130017.1941:(fields)
    #@+node:gcross.20090827130017.1942:__init__
    def __init__(self,slice_number,number_of_bins,filename):
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1942:__init__
    #@+node:gcross.20090827130017.1943:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        vpi.histograms.accumulate_neighbor_angular_separation_densities(
            angles,
            self.histogram
        )
    #@-node:gcross.20090827130017.1943:update
    #@-others
#@-node:gcross.20090827130017.1940:class NeighborAngularSeparationDensityHistogram
#@+node:gcross.20090827130017.1944:class AngularVelocityHistogram
class AngularVelocityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1945:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1945:__init__
    #@+node:gcross.20090827130017.1946:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        vpi.histograms.accumulate(first_derivatives,self.left,self.right,self.histogram)
    #@-node:gcross.20090827130017.1946:update
    #@-others
#@-node:gcross.20090827130017.1944:class AngularVelocityHistogram
#@+node:gcross.20090827130017.1947:class AngularVelocitySquaredHistogram
class AngularVelocitySquaredHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1948:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1948:__init__
    #@+node:gcross.20090827130017.1949:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        vpi.histograms.accumulate(first_derivatives,self.left,self.right,self.histogram)
    #@-node:gcross.20090827130017.1949:update
    #@-others
#@-node:gcross.20090827130017.1947:class AngularVelocitySquaredHistogram
#@+node:gcross.20090827130017.1950:class RotationQuadraticTermHistogram
class RotationQuadraticTermHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1951:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1951:__init__
    #@+node:gcross.20090827130017.1952:update
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
    #@-node:gcross.20090827130017.1952:update
    #@-others
#@-node:gcross.20090827130017.1950:class RotationQuadraticTermHistogram
#@+node:gcross.20090827130017.1953:class AngularVelocityAndRadiusHistogram
class AngularVelocityAndRadiusHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090827130017.1954:__init__
    def __init__(self,slice_number,maximum_angular_velocity,maximum_radius,number_of_bins,filename):
        self.maximum_angular_velocity = maximum_angular_velocity
        self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,)*2,dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090827130017.1954:__init__
    #@+node:gcross.20090827130017.1955:update
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
    #@-node:gcross.20090827130017.1955:update
    #@+node:gcross.20090827130017.1956:write_out_totals
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
    #@-node:gcross.20090827130017.1956:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1953:class AngularVelocityAndRadiusHistogram
#@-node:gcross.20090827130017.1910:Histograms
#@+node:gcross.20090827130017.1957:Energy estimates
#@+node:gcross.20090827130017.1958:class TotalEnergyEstimate
class TotalEnergyEstimate(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1960:update
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
    #@-node:gcross.20090827130017.1960:update
    #@-others
#@-node:gcross.20090827130017.1958:class TotalEnergyEstimate
#@+node:gcross.20090827130017.1962:Slice estimates
#@+node:gcross.20090827130017.1963:class SliceEnergyEstimate
class SliceEnergyEstimate(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.1964:__init__
    def __init__(self,slice_number,filename,label):
        SingleAverageValueEstimateAppendedToFile.__init__(self,filename,label)
        self.slice_number = slice_number
    #@-node:gcross.20090827130017.1964:__init__
    #@-others
#@-node:gcross.20090827130017.1963:class SliceEnergyEstimate
#@+node:gcross.20090827130017.1965:class EffectivePotentialSliceEnergyEstimate
class EffectivePotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1966:update
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
        self.add(sum(U))
    #@-node:gcross.20090827130017.1966:update
    #@-others
#@-node:gcross.20090827130017.1965:class EffectivePotentialSliceEnergyEstimate
#@+node:gcross.20090827130017.1967:class PhysicalPotentialSliceEnergyEstimate
class PhysicalPotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1968:update
    def update(self):
        self.add(sum(dot(self.system.x[self.slice_number]**2,self.system.harmonic_oscillator_coefficients))/2.0)
    #@-node:gcross.20090827130017.1968:update
    #@-others
#@-node:gcross.20090827130017.1967:class PhysicalPotentialSliceEnergyEstimate
#@+node:gcross.20090827130017.1969:class TotalPotentialSliceEnergyEstimate
class TotalPotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1970:update
    def update(self):
        self.add(sum(self.system.U[self.slice_number]))
    #@-node:gcross.20090827130017.1970:update
    #@-others
#@-node:gcross.20090827130017.1969:class TotalPotentialSliceEnergyEstimate
#@-node:gcross.20090827130017.1962:Slice estimates
#@+node:gcross.20090827130017.1971:Path estimates
#@+node:gcross.20090827130017.1972:class PathEnergyEstimates
class PathEnergyEstimates(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1973:(fields)
    estimates = 0
    #@-node:gcross.20090827130017.1973:(fields)
    #@+node:gcross.20090827130017.1974:__init__
    def __init__(self,filename):
        self.filename = filename
    #@-node:gcross.20090827130017.1974:__init__
    #@+node:gcross.20090827130017.1975:write_out_totals
    def write_out_totals(self,total_estimates):
        ensure_path_to_file_exists(self.filename)
        center_slice_number = self.system.center_slice_number
        with open(self.filename,"w") as f:
            for slice_number, estimate in enumerate(total_estimates):
                print >> f, center_slice_number-abs(center_slice_number-slice_number), estimate, slice_number
    #@-node:gcross.20090827130017.1975:write_out_totals
    #@-others
#@-node:gcross.20090827130017.1972:class PathEnergyEstimates
#@+node:gcross.20090827130017.1976:class EffectivePotentialPathEnergyEstimates
class EffectivePotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090827130017.1977:update
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
    #@-node:gcross.20090827130017.1977:update
    #@-others
#@-node:gcross.20090827130017.1976:class EffectivePotentialPathEnergyEstimates
#@+node:gcross.20090827130017.1978:class PhysicalPotentialPathEnergyEstimates
class PhysicalPotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090827130017.1979:update
    def update(self):
        self.estimates += sum(dot(self.system.x**2,self.system.harmonic_oscillator_coefficients),axis=-1)/2.0
    #@-node:gcross.20090827130017.1979:update
    #@-others
#@-node:gcross.20090827130017.1978:class PhysicalPotentialPathEnergyEstimates
#@+node:gcross.20090827130017.1980:class TotalPotentialPathEnergyEstimates
class TotalPotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090827130017.1981:update
    def update(self):
        self.estimates += sum(self.system.U,axis=-1)
    #@-node:gcross.20090827130017.1981:update
    #@-others
#@-node:gcross.20090827130017.1980:class TotalPotentialPathEnergyEstimates
#@-node:gcross.20090827130017.1971:Path estimates
#@-node:gcross.20090827130017.1957:Energy estimates
#@+node:gcross.20090827130017.1982:Position estimates
#@+node:gcross.20090827130017.1983:class AveragePositionEstimate
class AveragePositionEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    pass
#@-node:gcross.20090827130017.1983:class AveragePositionEstimate
#@+node:gcross.20090827130017.1984:class AverageAxialCoordinateEstimate
class AverageAxialCoordinateEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1985:__init__
    def __init__(self,axis,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        self.axis = axis
    #@-node:gcross.20090827130017.1985:__init__
    #@+node:gcross.20090827130017.1986:update
    def update(self):
        self.add(average(self.system.x[self.slice_number,:,self.axis]))
    #@-node:gcross.20090827130017.1986:update
    #@-others
#@-node:gcross.20090827130017.1984:class AverageAxialCoordinateEstimate
#@+node:gcross.20090827130017.1987:class AverageAxialDistanceEstimate
class AverageAxialDistanceEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1988:__init__
    def __init__(self,axis,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        self.axis = axis
    #@-node:gcross.20090827130017.1988:__init__
    #@+node:gcross.20090827130017.1989:update
    def update(self):
        self.add(average(abs(self.system.x[self.slice_number,:,self.axis])))
    #@-node:gcross.20090827130017.1989:update
    #@-others
#@-node:gcross.20090827130017.1987:class AverageAxialDistanceEstimate
#@+node:gcross.20090827130017.1990:class AverageRadiusEstimate
class AverageRadiusEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1991:update
    def update(self):
        self.add(vpi.observables.compute_radius_average(self.system.x[self.slice_number]))
    #@-node:gcross.20090827130017.1991:update
    #@-others
#@-node:gcross.20090827130017.1990:class AverageRadiusEstimate
#@+node:gcross.20090827130017.1992:class AveragePlaneRadiusEstimate
class AveragePlaneRadiusEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1993:__init__
    def __init__(self,plane_axis_1,plane_axis_2,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        assert plane_axis_1 >= 1
        assert plane_axis_2 >= 1
        assert not (plane_axis_1 == plane_axis_2)
        self.plane_axis_1 = plane_axis_1
        self.plane_axis_2 = plane_axis_2
    #@-node:gcross.20090827130017.1993:__init__
    #@+node:gcross.20090827130017.1994:update
    def update(self):
        self.add( vpi.observables.compute_plane_radius_average(self.system.x[self.slice_number],self.plane_axis_1,self.plane_axis_2))
    #@-node:gcross.20090827130017.1994:update
    #@-others
#@-node:gcross.20090827130017.1992:class AveragePlaneRadiusEstimate
#@+node:gcross.20090827130017.1995:class AverageRecipricalPlaneRadiusSquaredEstimate
class AverageRecipricalPlaneRadiusSquaredEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090827130017.1996:__init__
    def __init__(self,plane_axis_1,plane_axis_2,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        assert plane_axis_1 >= 1
        assert plane_axis_2 >= 1
        assert not (plane_axis_1 == plane_axis_2)
        self.plane_axis_1 = plane_axis_1
        self.plane_axis_2 = plane_axis_2
    #@-node:gcross.20090827130017.1996:__init__
    #@+node:gcross.20090827130017.1997:update
    def update(self):
        self.add(vpi.observables.compute_recip_plane_r_sq_average(self.system.x,self.plane_axis_1,self.plane_axis_2))
    #@-node:gcross.20090827130017.1997:update
    #@-others
#@-node:gcross.20090827130017.1995:class AverageRecipricalPlaneRadiusSquaredEstimate
#@-node:gcross.20090827130017.1982:Position estimates
#@+node:gcross.20090827130017.1998:Rotation related estimates
#@+node:gcross.20090827130017.1999:class AverageAngularVelocityEstimate
class AverageAngularVelocityEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.2000:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        self.add(average(first_derivatives))
    #@-node:gcross.20090827130017.2000:update
    #@-others
#@-node:gcross.20090827130017.1999:class AverageAngularVelocityEstimate
#@+node:gcross.20090827130017.2001:class AverageAngularVelocitySquaredEstimate
class AverageAngularVelocitySquaredEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.2002:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        self.add(average(first_derivatives))
    #@-node:gcross.20090827130017.2002:update
    #@-others
#@-node:gcross.20090827130017.2001:class AverageAngularVelocitySquaredEstimate
#@+node:gcross.20090827130017.2003:class AverageRotationQuadraticTermEstimate
class AverageAngularVelocitySquaredEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.2004:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpi.angular_momentum.compute_angular_derivatives(
            system.x[self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        self.add(average(first_derivatives))
    #@-node:gcross.20090827130017.2004:update
    #@-others
#@-node:gcross.20090827130017.2003:class AverageRotationQuadraticTermEstimate
#@+node:gcross.20090827130017.2008:class AverageAngularSeparationEstimate
class AverageAngularSeparationEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.2009:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        self.add(vpi.observables.compute_average_angular_separation(angles))
    #@-node:gcross.20090827130017.2009:update
    #@-others
#@-node:gcross.20090827130017.2008:class AverageAngularSeparationEstimate
#@+node:gcross.20090827130017.2010:class AverageNeighborAngularSeparationEstimate
class AverageNeighborAngularSeparationEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.2011:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        angles = arctan2(x[:,system.rotation_plane_axis_2],x[:,system.rotation_plane_axis_1])
        self.add(vpi.observables.compute_avg_neighbor_angular_sep(angles))
    #@-node:gcross.20090827130017.2011:update
    #@-others
#@-node:gcross.20090827130017.2010:class AverageNeighborAngularSeparationEstimate
#@+node:gcross.20090827130017.2012:class AverageRotationQuadraticTermEstimate
class AverageRotationQuadraticTermEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090827130017.2013:update
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
    #@-node:gcross.20090827130017.2013:update
    #@-others
#@-node:gcross.20090827130017.2012:class AverageRotationQuadraticTermEstimate
#@-node:gcross.20090827130017.1998:Rotation related estimates
#@-others
#@-node:gcross.20090827130017.1895:Observable classes
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
