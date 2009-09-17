<?xml version="1.0" encoding="utf-8"?>
<?xml-stylesheet ekr_test?>
<leo_file>
<leo_header file_format="2" tnodes="0" max_tnode_index="0" clone_windows="0"/>
<globals body_outline_ratio="0.5">
	<global_window_position top="22" left="1" height="564" width="1270"/>
	<global_log_window_position top="0" left="0" height="0" width="0"/>
</globals>
<preferences/>
<find_panel_settings/>
<vnodes>
<v t="gcross.20090728105516.1229" str_leo_pos="4"><vh>@chapters</vh></v>
<v t="gcross.20090728105516.1232"><vh>@ipy-startup</vh></v>
<v t="gcross.20090731093156.1320"><vh>nodeActions</vh>
<v t="gcross.20090731093156.1321"><vh>*.pdf</vh></v>
<v t="gcross.20090731093156.1322"><vh>*.dvi</vh></v>
<v t="gcross.20090731153600.1300"><vh>*.tex</vh></v>
<v t="gcross.20090731093156.1323"><vh>@files SConstruct</vh></v>
</v>
<v t="gcross.20090730140239.1300"><vh>Numerical analyses</vh>
<v t="gcross.20090728105516.1228"><vh>Proof that phase form is invalid</vh>
<v t="gcross.20090728110622.1243"><vh>via f</vh>
<v t="gcross.20090728110622.1241"><vh>f</vh></v>
<v t="gcross.20090728110622.1244"><vh>test run</vh></v>
</v>
<v t="gcross.20090728110622.1242"><vh>via fprim</vh>
<v t="gcross.20090728105516.1230"><vh>f'</vh></v>
<v t="gcross.20090728110622.1239"><vh>f'</vh></v>
<v t="gcross.20090728110622.1240"><vh>test run</vh></v>
</v>
</v>
<v t="gcross.20090729160450.1249"><vh>Check of backflow terms</vh>
<v t="gcross.20090729160450.1250"><vh>f</vh></v>
<v t="gcross.20090729160450.1252"><vh>g</vh></v>
<v t="gcross.20090729160450.1253"><vh>h</vh></v>
<v t="gcross.20090729160450.1254"><vh>pfun</vh></v>
<v t="gcross.20090729160450.1251"><vh>check</vh></v>
</v>
<v t="gcross.20090916084821.1352"><vh>Analysis of nodal surface</vh>
<v t="gcross.20090916085431.1357"><vh>Import needed modules</vh></v>
<v t="gcross.20090916084821.1353"><vh>Initialize system</vh></v>
<v t="gcross.20090916104053.1360"><vh>"Initialize system"</vh></v>
<v t="gcross.20090916085431.1356"><vh>Run optimization</vh></v>
</v>
<v t="gcross.20090916114839.1820"><vh>Analysis of linear correction thingie</vh>
<v t="gcross.20090916114839.1822"><vh>Import needed modules</vh></v>
<v t="gcross.20090916114839.1823"><vh>Initialize system</vh></v>
</v>
</v>
<v t="gcross.20090904201537.1554"><vh>@thin test_backflow.py</vh></v>
<v t="gcross.20090730140239.1301"><vh>LaTeX Notes</vh>
<v t="gcross.20090730140239.1305"><vh>@thin SConstruct</vh></v>
<v t="gcross.20090730140239.1304"><vh>@thin notes-prelude.in</vh></v>
<v t="gcross.20090730140239.3756"><vh>ExtraKineticTerms.pdf</vh>
<v t="gcross.20090730140239.1303"><vh>@thin ExtraKineticTerms.tex</vh></v>
</v>
<v t="gcross.20090817175122.1537"><vh>FixedPhase.pdf</vh>
<v t="gcross.20090817175122.1538"><vh>@thin FixedPhase.tex</vh></v>
</v>
<v t="gcross.20090730151355.1299"><vh>ExtraRotationTerms.pdf</vh>
<v t="gcross.20090730151355.1298"><vh>@thin ExtraRotationTerms.tex</vh></v>
</v>
<v t="gcross.20090730150937.1287"><vh>SPFixedPhaseApproximation.pdf</vh>
<v t="gcross.20090730140239.1306"><vh>@thin SPFixedPhaseApproximation.tex</vh></v>
</v>
<v t="gcross.20090730151355.1302"><vh>FancyFixedPhaseApproximation.pdf</vh>
<v t="gcross.20090730151355.1303"><vh>@thin FancyFixedPhaseApproximation.tex</vh></v>
</v>
<v t="gcross.20090730150937.1289"><vh>MagneticFieldAnalogy.pdf</vh>
<v t="gcross.20090730140239.1308"><vh>@thin MagneticFieldAnalogy.tex</vh></v>
</v>
<v t="gcross.20090803112359.1301"><vh>TermResumming.pdf</vh>
<v t="gcross.20090803112359.1303"><vh>@thin TermResumming.tex</vh></v>
</v>
<v t="gcross.20090827093130.1303"><vh>ImaginaryPart.pdf</vh>
<v t="gcross.20090827093130.1304"><vh>@thin ImaginaryPart.tex</vh></v>
</v>
<v t="gcross.20090904125219.1313"><vh>HardSphereTrialDerivatives.pdf</vh>
<v t="gcross.20090904125219.1316"><vh>@thin HardSphereTrialDerivatives.tex</vh></v>
</v>
<v t="gcross.20090914154930.2015"><vh>FancyPhaseGreensFunction.pdf</vh>
<v t="gcross.20090914154930.2017"><vh>@thin FancyPhaseGreensFunction.tex</vh></v>
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

import vpif
</t>
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
<t tx="gcross.20090729160450.1249"></t>
<t tx="gcross.20090729160450.1250">def f(r_i,r_j):
    return dot(
        (r_i-r_j)/norm(r_i-r_j),
        array((r_i[1],-r_i[0]))/norm(r_i)
    )
</t>
<t tx="gcross.20090729160450.1251">values = rand(2*2)

laplacian = sum(
    derivative(partial(pfun,h,n,values),values[n],dx=1e-6,n=2,order=13)
    for n in xrange(len(values))
    )

print laplacian
</t>
<t tx="gcross.20090729160450.1252">def g(coordinates):
    return sum(
        f(coordinates[i],coordinates[j])
        for i in xrange(len(coordinates))
        for j in xrange(len(coordinates))
        if i != j
    )
</t>
<t tx="gcross.20090729160450.1253">def h(values):
    return g(values.reshape((len(values)/2,2)))
</t>
<t tx="gcross.20090729160450.1254">def pfun(fn,n,values,new_value):
    old_value = values[n]
    values[n] = new_value
    result = fn(values)
    values[n] = old_value
    return result
</t>
<t tx="gcross.20090730140239.1300"></t>
<t tx="gcross.20090730140239.1301"></t>
<t tx="gcross.20090730140239.3756"></t>
<t tx="gcross.20090730150937.1287"></t>
<t tx="gcross.20090730150937.1289"></t>
<t tx="gcross.20090730151355.1299"></t>
<t tx="gcross.20090730151355.1302"></t>
<t tx="gcross.20090731093156.1320"></t>
<t tx="gcross.20090731093156.1321">import platform
if platform.system() == "Linux":
    g.os.system("kpdf " + pClicked.h + " &amp;")
else:
    g.os.system("open " + pClicked.h)
</t>
<t tx="gcross.20090731093156.1322">g.os.system("kdvi " + pClicked.h + " &amp;")
</t>
<t tx="gcross.20090731093156.1323">g.os.system("scons")

</t>
<t tx="gcross.20090731153600.1300">g.os.system("pdflatex " + " ".join(pClicked.h.split()[1:]) + " &amp;")
</t>
<t tx="gcross.20090803112359.1301"></t>
<t tx="gcross.20090817175122.1537"></t>
<t tx="gcross.20090827093130.1303"></t>
<t tx="gcross.20090904125219.1313"></t>
<t tx="gcross.20090914154930.2015"></t>
<t tx="gcross.20090916084821.1352"></t>
<t tx="gcross.20090916084821.1353">number_of_particles = 60
number_of_rotating_particles = 30

compute_amps_and_sum_syms_ampsq = vpif.angular_momentum.compute_amps_and_sum_syms_ampsq
def f(x):
    return compute_amps_and_sum_syms_ampsq(x.reshape(number_of_particles,2),number_of_rotating_particles,1,2)

compute_gradient_fancy_amplitude = vpif.angular_momentum.compute_gradient_fancy_amplitude
def fprime(x):
    return compute_gradient_fancy_amplitude(x.reshape(number_of_particles,2),number_of_rotating_particles,1,2).ravel()
</t>
<t tx="gcross.20090916085431.1356">x0 = rand(number_of_particles,2)
x0 /= sqrt(x0[:,0]**2 + x0[:,1]**2).reshape(number_of_particles,1) * 10

#x = fmin_bfgs(f,x0.ravel(),fprime)
#print f(x), norm(x0.ravel()-x)
#x,value,info = fmin_l_bfgs_b(f,x0.ravel(),None,(),True,[(x-0.1,x+0.1) for x in x0.ravel()],maxfun=100)
#print f(x), norm(x0.ravel()-x)
#x,_,value = fmin_tnc(f,x0.ravel(),None,(),True,[(x-1,x+1) for x in x0.ravel()],messages=0)
#print f(x), norm(x0.ravel()-x)

#x = fmin_cobyla(f,x0.ravel(),(lambda x: 5-norm(x-x0.ravel())))
#print f(x), norm(x-x0.ravel())
</t>
<t tx="gcross.20090916085431.1357">import vpif

from scipy.optimize import *
</t>
<t tx="gcross.20090916104053.1360">setup = """
import vpif

from numpy.random import rand

number_of_particles = 60
number_of_rotating_particles = 5

compute_amps_and_sum_syms_ampsq = vpif.angular_momentum.compute_amps_and_sum_syms_ampsq
def f(x):
    return compute_amps_and_sum_syms_ampsq(x.reshape(number_of_particles,2),number_of_rotating_particles,1,2)

compute_gradient_fancy_amplitude = vpif.angular_momentum.compute_gradient_fancy_amplitude
def fprime(x):
    return compute_gradient_fancy_amplitude(x.reshape(number_of_particles,2),number_of_rotating_particles,1,2).ravel()

x0 = rand(number_of_particles,2)
"""
</t>
<t tx="gcross.20090916114839.1820"></t>
<t tx="gcross.20090916114839.1822">from numpy.linalg import *

import vpif
</t>
<t tx="gcross.20090916114839.1823">number_of_particles = 10
number_of_rotating_particles = 5

x0 = rand(number_of_particles,2)

compute_amps_and_sum_syms = vpif.angular_momentum.compute_amps_and_sum_syms
def f(x):
    amplitude = compute_amps_and_sum_syms(x,number_of_rotating_particles,1,2)
    return array((real(amplitude),imag(amplitude)))

compute_gradient_fancy_amplitude = vpif.angular_momentum.compute_gradient_fancy_amplitude
def fprime(x):
    return compute_gradient_fancy_amplitude(x,number_of_rotating_particles)

print norm(lstsq(fprime(x0).transpose(),f(x0))[0])
</t>
</tnodes>
</leo_file>
