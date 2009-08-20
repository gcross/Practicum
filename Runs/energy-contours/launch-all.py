#@+leo-ver=4-thin
#@+node:gcross.20090818114910.1275:@thin launch-all.py
number_of_nodes = 8
wall_time = "0:02"

from numpy import arange

for frame_angular_momentum in [0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]:
    job_name = "vpi:%.2f" % frame_angular_momentum
    command = 'psub -b micphys -c atlas -i "srun -n 64 /g/g13/gcross/local/bin/python run.py {frame_angular_momentum} > data-3/run:{frame_angular_momentum}:.dat" -ln {number_of_nodes} -tM {wall_time} -r {job_name} -v -x'.format(**globals())
    print command
#@-node:gcross.20090818114910.1275:@thin launch-all.py
#@-leo
