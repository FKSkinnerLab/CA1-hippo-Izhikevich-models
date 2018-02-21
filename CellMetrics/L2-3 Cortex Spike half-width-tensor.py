# -*- coding: utf-8 -*-
#Created on Thu Nov 23 18:20:34 2017
#Modified on Feb 20, 2018
#Author: Amir
# Calculates the spike half-width of L2/3 Cortex cell model by iterating kHigh under control/4AP conditions


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


prefs.codegen.target = "numpy"
#start off each time
start_scope()

# Strongly adapting params (Izhikevitch)
Cm = 73 * pF # Membrane capacitance under control cond.
#Cm = 49 * pF # Membrane capacitance under 4AP cond.
vR = -60.6 * mV # Resting membrane potential
vPeak = 2.5 * mV # Spike cutoff
vT = -43.1 * mV # Instantaneous threshold membrane potential
 # Expected rheobase current under control cond.
I_rheoc = 43 * pA # Expected Rheobase current under cont. cond.
#I_rheoc = 26 * pA # Expected Rheobase current under 4AP cond.
Iappfixed = 400 * pA # Fixed Current Injection

# =============================================================================
# kLow = 0.3 * nS/mV
# kHigh = 14 * nS/mV
# =============================================================================

a = 0.1 /ms
d = 0.1 * pA
b = -0.2 * nS #Control
#b = -0.4 * nS #4AP
c = -67 * mV

#The values of the default current input - max and min. Sequence of events: Down-up-down
stepDown = 0.0 * pA
stepUp = 2 * pA

#Time breakdowns - the interval of applied current is between beforeStep and afterStep
beforeStep = 0.5 * second
afterStep = 1.5 * second

# kLow value iterations
klRes = 0.2 # 1.0E-4
klUp = 2
klLow = 0.0
klVals = [] #arange(klLow,klUp,klRes)

# kHigh iters
khRes = 1 # 0.5
khUp = 11
khLow = 1
khVals = arange(khLow,khUp,khRes)
InitShw = zeros((len(khVals),1)) # Initial Spike half-width
FinalShw = zeros((len(khVals),1)) # Final Spike half-width
ReqkH = zeros((2,len(khVals))) # Tensor of kHs
#Iterate over a and d
#for kli in range(len(klVals)):
for khi in range(len(khVals)):
    #Finding the parameter values in each iteration
    
    kLow = 0.6 * nS/mV #Control
    #kLow = 0.4 * nS/mV #4AP
    #klVals[kli] * nS/mV
    kHigh = khVals[khi] * nS/mV
    
    # Izhikevitch-Skinner cell equations
    IzhCelEq = '''
    # External current
# =============================================================================
#         Iapp = (t < beforeStep)*stepDown + (t >= beforeStep)*stepUp * i  + (t >= afterStep)*(-stepUp * i): amp
# =============================================================================
    k = (v<vT)*kLow+(v>=vT)*kHigh : siemens/volt
    du/dt = a*(b*(v-vR)-u) : amp
    dv/dt = (k*(v-vR)*(v-vT)+Iappfixed-u)/Cm : volt
    '''
    # Creating Cells
    Cells = NeuronGroup(1, model=IzhCelEq, reset = 'v = c; u += d' ,threshold='v>=vPeak', method='euler')
    #Resting potential used as base
    Cells.v = vR
    Cells.u = 0 * pA
    Cells_V = StateMonitor(Cells,'v',record=True) # Recording membrane potential for every neuron 
    Cells_time = SpikeMonitor(Cells)
    
    #Integration timecuts
    defaultclock.dt = 0.001*ms
    duration = 2000 * ms # Run time duration
    run(duration)
    PotenV = Cells_V.v[0] # Save Potential
    PotenT = Cells_V.t/ms # Save corresponding time steps for  potential
    SpikeT = Cells_time.t/ms
    # Half-width potential
    hwInitial = -20.3 * mV #((PotenV[ReturnIndex(PotenT,SpikeT[0])] + vT)/2)
#        hwFinal = -23.2 * mV # ((PotenV[ReturnIndex(PotenT,SpikeT[len(SpikeT)-1])] + vT)/2)
# =============================================================================
#         print('spike first %s' % SpikeT[0])
#         print('potent first %s' % PotenT[ReturnIndex(PotenV,hwInitial)])
#         print(ReturnIndex(Cells_V.t/ms,SpikeT[0]))
# =============================================================================
    PotenVLas = PotenV[(ReturnIndex(Cells_V.t/ms,SpikeT[len(SpikeT)-2])+1):(ReturnIndex(Cells_V.t/ms,SpikeT[len(SpikeT)-1])+1)]
    InitShw[khi,0] = (SpikeT[0] - PotenT[ReturnIndex(PotenV,hwInitial)])
    FinalShw[khi,0] = SpikeT[len(SpikeT)-1] - PotenT[ReturnIndex(PotenT,SpikeT[len(SpikeT)-2]) + ReturnIndex(PotenVLas,hwInitial)]#((vT+vPeak)/2))] # Final spike half-width
    print('Cm = %s, ' % Cm,'kLow = %s, ' % kLow,'kHigh = %s' % khVals[khi])
    print('Initial spike half-width = %s' %InitShw[khi,0])
    print('Final spike half-width = %s' %FinalShw[khi,0])
    figure(khi)
    plot(Cells_V.t/ms, Cells_V.v[0]/mV)
    xlabel("Time (ms)")
    ylabel("Potential (mV)")
    show()
    ReqkH[0,khi] = kHigh # Saving the kHighs as the first column of ReqkH tensor
    # Saving the spike half-widths which lies within 0.9 and 1.1 ms as 1 in the 2nd column of ReqkH and 0 otherwise
    ReqkH[1,khi] = (InitShw[khi,0]<1.1)*2 + (InitShw[khi,0]>1.1)*1 - (InitShw[khi,0]>0.9)*1 - (InitShw[khi,0]<0.9)*2 
print('Required kHigh values for Shw of 1 ms = %s' %ReqkH)

#Saving the Tensor of ReqkH
#save('Required-kHigh', ReqkH)
 
