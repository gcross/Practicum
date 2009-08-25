from math import pi
import os

for total_number_of_particles in map(int,os.listdir("data")):
    for frame_angular_velocity in os.listdir("data/{total_number_of_particles}".format(**vars())):
        with open("data/{total_number_of_particles}/{frame_angular_velocity}/average-radii".format(**vars()),"w") as out:
            number_of_bins = 51
            maximum_radius = 2.5
            bin_width = maximum_radius/number_of_bins
            half_bin_width = bin_width
            start = half_bin_width
            for n in xrange(0,total_number_of_particles+1):
                average = 0
                counts = 0
                current = start
                try:
                    with open("data/{total_number_of_particles}/{frame_angular_velocity}/{n}/center/radial-density".format(**vars()),"r") as f:
                        for line in f:
                            bin, count = map(float,line.split())
                            normalization = (4.0/3.0)*pi*((current+half_bin_width)**3-(current-half_bin_width)**3)
                            count *= normalization

                            average += bin*count
                            counts += count

                            current += bin_width
                    print >> out, n, average/counts
                except:
                    pass
