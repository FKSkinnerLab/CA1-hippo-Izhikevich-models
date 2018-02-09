# -*- coding: utf-8 -*-
#Created on Nov 9, 2017
#Modified on Feb 8, 2018
#Author: Amir
#Generates the single cell PV+ model of Katie Work in 2013 

from brian2 import *
from pylab import *

# Calling/recalling the time index of called potential/time
def ReturnIndex(original,ref):
    hh = 0
    for h in range(len(original)):
        if (original[h] == ref):
            return h
        elif (original[h] <= ref):
            hh = h + 1
        else: return hh-1

#No compiling - remove if needed
prefs.codegen.target = "numpy"
#Don't include previous Brian objects
start_scope()

# Strongly adapting params (Izhikevitch)
Cm = 90 * pF # Membrane capacitance
vR = -60.6 * mV # Resting membrane potential
vPeak = 2.5 * mV # Spike cutoff
vT = -43.1 * mV # Instantaneous membrane potential
Ishift = 131.0 * pA # Rheobase current shift
Iother = 0.0 * pA # Defined at zero since input is not noisy
IappFixed = 400 * pA
N = 60

kLow = 1.7 * nS/mV
kHigh = 14 * nS/mV

a = 0.1 /ms
d = 0.1 * pA
b = -0.1 * nS
c = -67 * mV

# Varying/Fixing the current input
Ivaried = arange(480,496,2)
IvariedF = 260 * pA #Ivaried[ii] *pA
    

#Izhikevitch-Skinner cell equations
cellEqs = '''
#Excitatory term
k = (v<=vT)*kLow+(v>vT)*kHigh : siemens/volt
du/dt = (a*(b*(v-vR)-u)) : amp
dv/dt = (k*(v-vR)*(v-vT)+IvariedF-u)/Cm : volt
'''
#Create cells
CELLS = NeuronGroup(1, model=cellEqs, threshold="v >= vPeak",  reset ="v = c; u += d" , method="euler")
#set initial conditions
#Resting potential used as base
CELLS.v = vR

#Record data
Cells_I = StateMonitor(CELLS,'u',record=0)
Cells_V = StateMonitor(CELLS,'v', record=0)
Cells_time = SpikeMonitor(CELLS) #Don't record to speed up runtime

 #Integration timecuts
defaultclock.dt = 0.001*ms

 #Set duration and run simulation
duration = 100 * ms
run(duration)

PotenV = Cells_V.v[0] # Save Potential
PotenT = Cells_V.t/ms # Save corresponding time steps for  potential
print('time steps = %s' % len(PotenT)) # Print number of time steps
print('current input = %s pA' % (IvariedF/pA))
SpikeT = Cells_time.t/ms

# =============================================================================
# # For unknown vT and vPeak, calculating the membrane potential of spike half-width
# HVinitial = ((PotenV[ReturnIndex(PotenT,SpikeT[0])] + vT)/2)
# HVfinal = ((PotenV[ReturnIndex(PotenT,SpikeT[len(SpikeT)-1])] + vT)/2)
# =============================================================================

PotenVLas = PotenV[ReturnIndex(PotenT,SpikeT[len(SpikeT)-2])+1:ReturnIndex(PotenT,SpikeT[len(SpikeT)-1])+1] # Potentials of final spike width
SWI = SpikeT[0] - PotenT[ReturnIndex(PotenV,vT)] # Initial spike width
print('initial spike width = %s ms' % SWI)
SWF = SpikeT[len(SpikeT)-1] - PotenT[ReturnIndex(PotenT,SpikeT[len(SpikeT)-2]) + ReturnIndex(PotenVLas,vT)] # Final spike width
print('final spike width = %s ms' % SWF)

# =============================================================================
# Finding the current input of called spike width
#     if (SWI < 0.84):
#         if (SWI > 0.82):
#             print(IvariedF)
# =============================================================================

 #Ploting the mv vs time for the set values
figure(1)
plot(Cells_V.t/ms, Cells_V.v[0]/mV)
xlabel ('Time (ms)')
ylabel('Potential (mV)')
