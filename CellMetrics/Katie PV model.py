# -*- coding: utf-8 -*-
#Created on Nov 9, 2017
#Author: Amir
#Generates a single cell PV+ model from Katie Works 2013 

from brian2 import *
from pylab import *

#No compiling - remove if needed
prefs.codegen.target = "numpy"
#Don't include previous Brian objects
start_scope()

# Strongly adapting params (Izhikevitch)
Cm = 90 * pF # Membrane capacitance
vR = -60.6 * mV # Resting membrane potential
vPeak = 2.5 * mV # Spike cutoff
vT = -43.1 * mV # Instantaneous membrane potential
Ishift = 0.0 * pA # Rheobase current shift
Iother = 0.0 * pA # Defined at zero since input is not noisy
IappFixed = 260 * pA

kLow = 1.7 * nS/mV
kHigh = 14 * nS/mV

a = 0.1 /ms
d = 0.1 * pA
b = -0.1 * nS
c = -67 * mV

                #Izhikevitch-Skinner cell equations
cellEqs = '''
                #Excitatory term
# =============================================================================
#                 Iapp = (t < beforeStep)*stepDown + (t >= beforeStep)*stepUp * i  + (t >= afterStep)*(-stepUp * i): amp
# =============================================================================
                k = (v<vT)*kLow+(v>=vT)*kHigh : siemens/volt
                du/dt = (a*(b*(v-vR)-u)) : amp
                dv/dt = ((k*(v-vR)*(v-vT)+IappFixed-u)/Cm) : volt 
                '''

                #Create cells
CELLS = NeuronGroup(1, model=cellEqs, reset ="v = c; u += d" , threshold="v>= vPeak", method="euler")
#set initial conditions
#Resting potential used as base
CELLS.v = vR


                #Record data
# =============================================================================
# Cells_I = StateMonitor(CELLS,'u',record=0)
# =============================================================================
Cells_V = StateMonitor(CELLS,'v', record= 0)
Cells_time = SpikeMonitor(CELLS) #Don't record to speed up runtime


                #Set duration and run simulation

duration = 100 * ms
run(duration)
 
print(Cells_V.v[0]/mV)
print(Cells_V.t/ms)         
                #Ploting the mv vs time for the set values

plot(Cells_V.t/ms, Cells_V.v[0]/mV)
