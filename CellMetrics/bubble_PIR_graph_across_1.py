# -*- coding: utf-8 -*-
#Created on June 23, 2017
#Modified on June 26, 2017
#Author: Anton Lunyov

from brian2 import *
from pylab import *

#No compiling to speed up process
prefs.codegen.target = "numpy"
#Don't include previous Brian objects
start_scope()

# Strongly adapting params (Izhikevitch)
Cm = 115.0 * pF # Membrane capacitance
vR = -61.8 * mV # Resting membrane potential
vPeak = 22.6 * mV # Spike cutoff
vT = -57.0 * mV # Instantaneous membrane potential
Ishift = 0.0 * pA # Rheobase current shift
Iother = 0.0 * pA # Defined at zero since input is not noisy

kLow = 0.1 * nS/mV
kHigh = 3.3  * nS/mV

#Izhikevitch params
a = 0.0012 /ms
d = 10.0 * pA
b = 3.0 * nS
c = -65.8 * mV

#Single neuron simulation - 20 neurons with varying current steps
N = 20
freq_initial = zeros(N)
freq_final = zeros(N)
Erev = -15 * mV #Reversal potential of the synapse, -15?

#The values of the default current input - max and min. Down-up-down
stepDown = 0.0 * pA
stepUp = -2.0 * pA

#Time breakdowns - the interval of applied current is between beforeStep and afterStep
beforeStep = 0.5 * second
afterStep = 1.5 * second

# a value iterations
aRes = 1.2E-4
aUp = 2.4E-3
aLow = 0.0E-3
aVals = arange(aLow,aUp,aRes)

# klow iterations
kRes = 0.05
kLowUp = 1.0
kLowLower = 0.0
kVals = arange(kLowLower,kLowUp,kRes)

# b iters
bRes = 0.3
bUp = 12.0
bLow = 0.0
bVals = arange(bLow,bUp,bRes)

# d iters
dRes = 0.5
dUp = 20
dLow = 0
dVals = arange(dLow,dUp,dRes)


# Array to store the currents at which a post-inhibitory spike is seen
transitionCurrents = zeros(len(aVals))

for y in range(len(aVals)):
    a = aVals[y] /ms

    #Izhikevitch-Skinner cell equations
    cellEqs = """
    #Excitatory term
    Iapp = (t < beforeStep)*stepDown + (t >= beforeStep)*stepUp * i  + (t >= afterStep)*(-stepUp * i): amp

    k = (v<vT)*kLow+(v>=vT)*kHigh : siemens/volt
    du/dt = a*(b*(v-vR)-u) : amp
    dv/dt = (k*(v-vR)*(v-vT)+Iapp-u)/Cm : volt 
    """

    #Create cells
    CELLS = NeuronGroup(N, model=cellEqs, reset ="v = c; u += d" , threshold="v>=vPeak", method="euler")

    #set initial conditions
    #Resting potential used as base
    CELLS.v = vR


    #Record plotted data
    Cells_V = StateMonitor(CELLS,'v',record=True)
    Spiketimes = SpikeMonitor(CELLS)

    defaultclock.dt = 0.1*ms

    #Set duration and run simulation
    duration = 2 * second

    run(duration)

    #These neurons are spiking before the membrane potential is decreased, and are thus ineligible for PIR analysis (in my limited expertise)
    falseSpikes = len(Spiketimes.t[Spiketimes.t < 1.5*second]/second)

    #If was not spiking before,
    if (falseSpikes == 0):
        spikes = zeros(N)
        for h in range(N):
            #Record the number of spikes after the hyperpolarization is over
            foo = Spiketimes.t[Spiketimes.i == h]
            foo = foo[foo >= 1.5*second]/second
            spikes[h] = len(foo)

        #Check where the transition occurs, put one zero so array is not empty (checks are easier)
        transitions = zeros(1)

	#If at zero point, save that as transition current
        if (spikes[0] > 1):
            transitions = append(transitions,0.0)

        #Otherwise find the exact point where the transition occurs. Choose latter as of the two as the point
        for k in range(1,len(spikes)):
            if ((int(spikes[k]) > 0) & (int(spikes[k - 1]) == 0)):
                transitions = append(transitions,k*stepUp/pA)
        
	#Eliminate the zero inserted in beginning, and record number of transitions
        transitionNum = len(transitions) - 1
		
	#Invalid TC = +1 pA
        #If more than one transitions occured, the effect is invalid - model breaks down
        if (transitionNum > 1):
            transitionCurrents[y] = +1
        #No transition occured - current required is not enough, thus invalid
        elif (transitionNum == 0):
            transitionCurrents[y] = +1
        elif (transitionNum == 1):
            transitionCurrents[y] = transitions[1]
    #If cell firing before, invalid
    else:
        transitionCurrents[y] = +1


print("Graphing")

#Plot the aggregate data
figure(1)    
subplot(211)
plot(aVals,transitionCurrents,'ob')
xlabel("a")
ylabel("Transition current (pA)")


#Plot potential for one cell
cell = 10
subplot(212)
plot(Cells_V.t/ms,Cells_V.v[cell]/mV) 
xlabel("Time (ms)")
ylabel("Cell 10 - PYR Membrane Potential (mV)")

show()
