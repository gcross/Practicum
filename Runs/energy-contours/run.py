#@+leo-ver=4-thin
#@+node:gcross.20090818114910.1235:@thin run.py
import sys
sys.path.append("lib")

from probe import probe, my_rank

n_particles = int(sys.argv[1])
frame_angular_velocity = float(sys.argv[2])

if my_rank == 0:
    for N_rotating_particles in xrange(n_particles+1):
        energy, denergy = probe(n_particles,frame_angular_velocity,N_rotating_particles)
        print N_rotating_particles, energy, denergy
else:
    for N_rotating_particles in xrange(n_particles+1):
        probe(n_particles,frame_angular_velocity,N_rotating_particles)
#@-node:gcross.20090818114910.1235:@thin run.py
#@-leo
