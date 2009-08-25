from math import pi
from contextlib import nested
from itertools import izip
import os

for total_number_of_particles in map(int,os.listdir("data")):
    for frame_angular_velocity in os.listdir("data/{total_number_of_particles}".format(**vars())):
        with nested(
            open("data/{total_number_of_particles}/{frame_angular_velocity}/kinetic-energy".format(**vars()),"w"),
            open("data/{total_number_of_particles}/{frame_angular_velocity}/center/total-potential".format(**vars()),"r"),
            open("data/{total_number_of_particles}/{frame_angular_velocity}/total-energy".format(**vars()),"r"),
        ) as (k_file,p_file,t_file):
            for (p_line, t_line) in izip(p_file,t_file):
                p_bin, p_energy = p_line.split()
                t_bin, t_energy, t_diff = t_line.split()
                assert (p_bin == t_bin)
                k_bin = p_bin
                k_energy = float(t_energy) - float(p_energy)
                print >> k_file, k_bin, k_energy
