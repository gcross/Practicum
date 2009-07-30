#@+leo-ver=4-thin
#@+node:gcross.20090730140239.1305:@thin SConstruct
#@@language Python

import os

#### Some configurations.

DEFAULT_TARGET = 'dvi'

#### Initialization

env = Environment(
    tools = ['tex','pdftex'],
    ENV = os.environ
    )

#### The actual builds.

## LaTeX DVI build:
for name in [
    "ExtraKineticTerms",
    "ExtraRotationTerms",
    "FancyFixedPhaseApproximation",
    "SPFixedPhaseApproximation",
    "MagneticFieldAnalogy",
  ]:
    dviOutput = env.DVI(source=name+'.tex', target=name+'.dvi')
    env.Alias('dvi', name+'.dvi')

    pdfOutput = env.PDF(source=name+'.tex', target=name+'.pdf')
    env.Alias('pdf', name+'.pdf')

Default(DEFAULT_TARGET)
#@-node:gcross.20090730140239.1305:@thin SConstruct
#@-leo
