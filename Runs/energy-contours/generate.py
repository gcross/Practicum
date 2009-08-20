for n_particles in [5,6]:
    for frame_angular_velocity in [0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]:
        print "echo Now writing: data-{n_particles}/run:{frame_angular_velocity}:.dat".format(**vars())
        print "srun -n 128 /g/g13/gcross/local/bin/python run.py {n_particles} {frame_angular_velocity} | tee data-{n_particles}/run:{frame_angular_velocity}:.dat".format(**vars())
