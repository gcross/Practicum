#@+leo-ver=4-thin
#@+node:gcross.20090818114910.1235:@thin run.py
import sys
sys.path.append("lib")

from probe import probe

frame_angular_velocity = float(sys.argv[1])

for N_rotating_particles in xrange(5):
    probe("run",frame_angular_velocity,N_rotating_particles)
#@-node:gcross.20090818114910.1235:@thin run.py
#@-leo
