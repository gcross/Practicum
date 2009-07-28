<?xml version="1.0" encoding="utf-8"?>
<?xml-stylesheet ekr_test?>
<leo_file>
<leo_header file_format="2" tnodes="0" max_tnode_index="0" clone_windows="0"/>
<globals body_outline_ratio="0.5">
	<global_window_position top="10" left="470" height="600" width="800"/>
	<global_log_window_position top="0" left="0" height="0" width="0"/>
</globals>
<preferences/>
<find_panel_settings/>
<vnodes>
<v t="gcross.20090728105516.1229" str_leo_pos="3,0,1"><vh>@chapters</vh></v>
<v t="gcross.20090728105516.1232"><vh>@ipy-startup</vh></v>
<v t="gcross.20090728105516.1233" a="E"><vh>@settings</vh>
<v t="gcross.20090728105516.1234"><vh>@string ipython_argv = ipython -pylab</vh></v>
</v>
<v t="gcross.20090728105516.1228" a="E"><vh>Proof that phase form is valid</vh>
<v t="gcross.20090728110622.1243" a="E"><vh>via f</vh>
<v t="gcross.20090728110622.1241"><vh>f</vh></v>
<v t="gcross.20090728110622.1244"><vh>test run</vh></v>
</v>
<v t="gcross.20090728110622.1242" a="E"><vh>via fprim</vh>
<v t="gcross.20090728105516.1230"><vh>f'</vh></v>
<v t="gcross.20090728110622.1239"><vh>f'</vh></v>
<v t="gcross.20090728110622.1240"><vh>test run</vh></v>
</v>
</v>
</vnodes>
<tnodes>
<t tx="gcross.20090728105516.1228"></t>
<t tx="gcross.20090728105516.1229"></t>
<t tx="gcross.20090728105516.1230">def fprime(n,angles):
    numerator = sum(cos(angles-angles[n,...]),axis=0)
    angles_1 = angles.reshape((1,) + angles.shape)
    angles_2 = angles.reshape(angles.shape[:1] + (1,) + angles.shape[1:])
    denominator = sum(sum(cos(angles_1-angles_2),axis=0),axis=0)
    return numerator/denominator
</t>
<t tx="gcross.20090728105516.1232">from numpy import *
from scipy.misc import *
from functools import *
</t>
<t tx="gcross.20090728105516.1233"></t>
<t tx="gcross.20090728105516.1234"></t>
<t tx="gcross.20090728110622.1239">def fprime(n,old_angles,replacement_angles):
    angles = old_angles.copy()
    angles[n] = replacement_angles
    return (
      sum(
        (cos(angles[i]-angles[n])
            for i in xrange(angles.shape[0])
        ),axis=0
        )
      /
      sum(
        (cos(angles[i]-angles[j])
            for i in xrange(angles.shape[0])
            for j in xrange(angles.shape[0])
        ),axis=0
        )
    )</t>
<t tx="gcross.20090728110622.1240">angles = rand(10)*2*pi

laplacian = 0

for n in xrange(angles.shape[0]):
    fpn = partial(fprime,n,angles)
    laplacian += derivative(fpn,angles[n],dx=1e-10)

print laplacian
</t>
<t tx="gcross.20090728110622.1241">def f(n,old_angles,replacement_angles):
    angles = old_angles.copy()
    angles[n] = replacement_angles
    return sum(sin(angles),axis=0)/sum(cos(angles),axis=0)
</t>
<t tx="gcross.20090728110622.1242"></t>
<t tx="gcross.20090728110622.1243"></t>
<t tx="gcross.20090728110622.1244">angles = rand(10)*2*pi

laplacian = 0

for n in xrange(angles.shape[0]):
    fn = partial(f,n,angles)
    laplacian += derivative(fpn,angles[n],dx=1e-6,n=2,order=13)

print laplacian
</t>
</tnodes>
</leo_file>
