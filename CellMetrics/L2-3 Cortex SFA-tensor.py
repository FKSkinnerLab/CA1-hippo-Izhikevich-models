# -*- coding: utf-8 -*-
# Created on Sat Nov 18 22:26:51 2017
# Modified on Thr Feb 20 19:00:00 2018
# author: Amir
# Itiration over a and d to achieve minimal SFA

from brian2 import *
from pylab import *

prefs.codegen.target = "numpy"
#start off each time
start_scope()

# Strongly adapting params (Izhikevitch)

Cm = 73 * pF # Membrane capacitance under control cond.
# Cm = 49 * pF # Membrane capacitance under 4AP cond.
vR = -60.6 * mV # Resting membrane potential
vPeak = 2.5 * mV # Spike cutoff
vT = -43.1 * mV # Instantaneous threshold membrane potential
 
I_rheoc = 43 * pA # Expected rheobase current under control cond.
#I_rheo4ap = 26 * pA  # Expected rheobase current under 4AP cond.

# Fixed for Control
kLow = 0.6 * nS/mV
kHigh = 1 * nS/mV

# Fixed for 4AP
#kLow = 0.4 * nS/mV
#kHigh = 0.65 * nS/mV

# Rest of parameters for hippocampul cell model
# =============================================================================
# a = 0.1 /ms
# d = 0.1 * pA
# =============================================================================
c = -67 * mV
b = -0.2 * nS #control
#b = -0.4 * nS #4AP
# Single neuron simulation - Number of current input
N = 60

#The values of the default current input - max and min. Sequence of events: Down-up-down
stepDown = 0.0 * pA
stepUp = 2 * pA

#Time breakdowns - the interval of applied current is between beforeStep and afterStep
beforeStep = 0.5 * second
afterStep = 1.5 * second

# a value iterations
aRes = 0.01 # 1.0E-4
aUp = 0.1
aLow = 0.0
aVals = arange(aLow,aUp,aRes)

# d iters
dRes = 0.25 # 0.5
dUp = 2.5
dLow = 0
dVals = arange(dLow,dUp,dRes)
#print(dVals)
#Contains ISI data
freq_initial = zeros(N)
freq_final = zeros(N)

#Tensors to store all data obtained - 10x10x10x10
adaptation = zeros((len(bVals)+1,len(kVals)+1))
adaptationRatio = zeros((len(bVals)+1,len(kVals)+1))
Irheobk = zeros((len(bVals),len(kVals)))

#Iterate over b and kLow
for ai in range(len(aVals)):
    for di in range(len(dVals)):
        count = 0 # counting the rheobase current
        # Finding the parameter values in each iteration
        a= aVals[ai] /ms
        d = dVals[di] * pA
        # Izhikevitch-Skinner cell equations
        IzhCelEq = '''
        # External current
# =============================================================================
#         Iapp = (t < beforeStep)*stepDown + (t >= beforeStep)*stepUp * i  + (t >= afterStep)*(-stepUp * i): amp
# =============================================================================
        Iapp = stepUp * i : amp
        k = (v<vT)*kLow+(v>=vT)*kHigh : siemens/volt
        du/dt = a*(b*(v-vR)-u) : amp
        dv/dt = (k*(v-vR)*(v-vT)+Iapp-u)/Cm : volt
        '''
        # Creating Cells
        Cells = NeuronGroup(N, model=IzhCelEq, reset = 'v = c; u += d' ,threshold='v>=vPeak', method='euler')
        #Resting potential used as base
        Cells.v = vR
        Cells.u = 0 * pA
        Cells_V = StateMonitor(Cells,'v',record=True) # Recording membrane potential for every neuron 
        SpikeT = SpikeMonitor(Cells)
        
        #Integration timecuts
        defaultclock.dt = 0.01 * ms
        duration = 2 * second # Run time duration
        run(duration)
        for h in range(N):
            # Retrieve the number of spikes for every neuron (and thus every current step)
            spikes = SpikeT.t[SpikeT.i == h]
            spikeNum = len(spikes)
            #Record ISI based on conditions
            if (spikeNum == 0):
                freq_initial[h] = 0
                freq_final[h] = 0
            elif (spikeNum == 1):
                freq_initial[h] = 1
                freq_final[h] = 1
            else:
                # First and last ISI
                freq_initial[h] = 1/(spikes[1] - spikes[0])                            
                freq_final[h] = 1/(spikes[spikeNum - 1] - spikes[spikeNum - 2])
        # Re-create the implicit x array of currents to fit to
        x = range(0,int(stepUp*N/pA),int(stepUp/pA))
        #Fit initial and final adaptation data
        fitInitial = polyfit(x,freq_initial,1)
        slopeInitial = fitInitial[0]
        fitFinal = polyfit(x,freq_final,1)
        slopeFinal = fitFinal[0]
        #Calculate relevant params
        SFA = slopeInitial - slopeFinal
        SFA_rat = (slopeInitial/slopeFinal)
        adaptation[ai+1,di+1] = SFA
        adaptationRatio[ai+1,di+1] = SFA_rat
        
        # Assigning the 1st row and column for the a and d values
        adaptation[0,di] = d
        adaptation[ai,0] = a
        adaptationRatio[0,di] = d
        adaptationRatio[ai,0] = a
        
            # Real-time printing the data with the SFA results for visual tracking    
        print('ai= %s' % a, 'di= %s' % d)
        print('SFA = %s' %SFA, 'SFA ratio = %s' %SFA_rat)

# Viewing the result for adaptation and its ratio
print('adaptation tensor presented as "a" rows by "d" columns')
print (adaptation)
print(adaptationRatio)

#Saving the adaptation tensors in npy files
save('adaptation_control', adaptation)
save('adaptationRatio_control', adaptationRatio)
