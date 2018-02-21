# -*- coding: utf-8 -*-
# Created on Sat Nov 18 22:26:51 2017
# Modified on Thr Feb 20 19:00:00 2018
# author: Amir
# Itiration over b and kLow to achieve the desired Irheo & F-I curve

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

# Assumed kHigh
#kLow = 0.3 * nS/mV
kHigh = 14 * nS/mV

# Rest of parameters for hippocampul cell model
a = 0.1 /ms
d = 0.1 * pA
c = -67 * mV
#b = -0.1 * nS

# Single neuron simulation - Number of current inputs
N = 60

#The values of the default current input - max and min. Sequence of events: Down-up-down
stepDown = 0.0 * pA
stepUp = 2 * pA

#Time breakdowns - the interval of applied current is between beforeStep and afterStep
beforeStep = 0.5 * second
afterStep = 1.5 * second

# klow iterations
kRes = 0.1 # 0.01
kLowUp = 1
kLowLower = 0.00
kVals = arange(kLowLower,kLowUp,kRes)

# b iters
bRes = 0.1 #0.2
bUp = 0.0
bLow = -1.0
bVals = arange(bLow,bUp,bRes)

#Contains ISI data
freq_initial = zeros(N)
freq_final = zeros(N)

#Tensor to store the obtained rheobase currents - 10x10
Irheobk = zeros((len(bVals),len(kVals)))

#Iterate over b and kLow
for bi in range(len(bVals)):
    for ki in range(len(kVals)):
        count = 0 # counting the rheobase current
        # Finding the parameter values in each iteration
        b = bVals[bi] *nS
        kLow = kVals[ki] * nS/mV
        # Izhikevitch-Skinner cell equations
        IzhCelEq = '''
        # External current
        Iapp = (t < beforeStep)*stepDown + (t >= beforeStep)*stepUp * i  + (t >= afterStep)*(-stepUp * i): amp
        #Iapp = stepUp * i : amp
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
                
        defaultclock.dt = 0.01 * ms #Integration timecuts
        duration = 1.5 * second # Run time duration
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
                x = arange(0,int(stepUp*N/pA),int(stepUp/pA))
            
            # Counting the number of input currents with zero frequency so as to calculate the Irheo
            count += (freq_initial[h] < 1) * 1 #+ (freq_initial[h] > 15) * 0   
# =============================================================================
#         fitInitialRheo = polyfit(x[count-1:count+4],freq_initial[count-1:count+4],1)
#         IrheobkFit [bi,ki] = ((-1)*fitInitialRheo[1])/fitInitialRheo[0]
# =============================================================================
               
        Irheobk [bi,ki] = x[count-1] # Savng the rheobase currents for indiv. comb. of b & kLow
# =============================================================================
#         
#       # Real-time visual tracking of the FI curve results
        print('Cm = %s, ' %Cm,'kLow = %s, ' %kVals[ki], 'b = %s' %bVals[bi])
        print('Input Current Index = %s' %count)
        print('Rheobase I = %s' %(Irheobk [bi,ki]))
        figure(ki)
        plot(x,freq_initial,'g')
        xlabel("Input Current (pA)")
        ylabel("Frequency (Hz)")
        show()

# Saving the tensor of rheobase current
save('Rheo-I-b-by-kLow', Irheobk)

