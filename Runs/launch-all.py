#@+leo-ver=4-thin
#@+node:gcross.20090818114910.1275:@thin launch-all.py
number_of_nodes = 8
wall_time = "0:02"

from numpy import arange

l = map(lambda x: x/10.0, xrange(1,10)) + [2,3,4,5,6,7,8,9,10]

for frame_angular_momentum in l+[0]+map(lambda x: -x,l):
    job_name = "run-ang-%.2f" % frame_angular_momentum
    command = 'psub -b micphys -c zeus -i "srun -n 64 /g/g13/gcross/local/bin/python run.py {frame_angular_momentum}" -ln {number_of_nodes} -tM {wall_time} -r {job_name} -v -x'.format(**globals())
    print command
#@-node:gcross.20090818114910.1275:@thin launch-all.py
#@-leo
