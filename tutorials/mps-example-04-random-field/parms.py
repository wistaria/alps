# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2016 by Michele Dolfi <dolfim@phys.ethz.ch>
# 
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot

basename = 'parm_random'

#prepare the input parameters
parms = []
for seed in [1,10,20]:
    parms.append( { 
        'LATTICE'                               : "open chain lattice",
        'L'                                     : 10,
        'MODEL'                                 : "spin",
        'J'                                     : 1,
        'Jz'                                    : 1,
        'h'                                     : 'gaussian_random(2,1)',
        'SWEEPS'                                : 5,
        'DISORDERSEED'                          : seed,
        'NUMBER_EIGENVALUES'                    : 1,
        'MAXSTATES'                             : 100,
        'MEASURE_LOCAL[Local mag]'              : 'Sz'
       } )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles(basename,parms)
res = pyalps.runApplication('mps_optim',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix=basename), ['Local mag'])

for d in pyalps.flatten(data):
    d.y = d.y[0]
    d.props['line'] = '-o'
    d.props['label'] = 'seed={}'.format(d.props['DISORDERSEED'])


plt.figure()
pyalps.plot.plot(data)
plt.legend()
plt.ylabel('local magnetization')
plt.xlabel('site')

plt.show()
