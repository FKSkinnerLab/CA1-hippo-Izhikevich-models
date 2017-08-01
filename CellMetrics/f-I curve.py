# -*- coding: utf-8 -*-
#Created on June 2, 2017
#Modified on June 2, 2017
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

a = 0.0012 /ms
d = 10.0 * pA
b = 3.0 * nS
c = -65.8 * mV

#Single neuron simulation - 20 neurons with varying current steps
N = 20
#Arrays for storing the initial and final ISI
freq_initial = zeros(N)
freq_final = zeros(N)

#The values of the default current input - max and min. Down-up-down
stepDown = 0.0 * pA
stepUp = 5.0 * pA

#Time breakdowns - the interval of applied current is between beforeStep and afterStep
beforeStep = 0.5 * second
afterStep = 1.5 * second

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

#For every neuron
for h in range(N):
    #Retrieve the number of spikes
    spikes = Spiketimes.t[Spiketimes.i == h]
    spikeNum = len(spikes)

    #Record ISI based on conditions
    if (spikeNum == 0):
        freq_initial[h] = 0
        freq_final[h] = 0
    elif (spikeNum == 1):
        freq_initial[h] = 1
        freq_final[h] = 1
    else:
        #First and last ISI
        freq_initial[h] = 1/(spikes[1] - spikes[0])
        freq_final[h] = 1/(spikes[spikeNum - 1] - spikes[spikeNum - 2])

#Do a linear fit of the data
x = arange(0,stepUp*N/pA,stepUp/pA)
fitInitial = polyfit(x,freq_initial,1)
fitFinal = polyfit(x,freq_final,1)

fitInitPlot = (fitInitial[0] * array(x)) + fitInitial[1]
fitFinalPlot = (fitFinal[0] * array(x)) + fitFinal[1]


#Plot relevant data
figure(1)    
plot(x,freq_initial,'ob')
plot(x,fitInitPlot,'b')
plot(x,freq_final,'og')
plot(x,fitFinalPlot,'g')
xlabel("I")
ylabel("f")
show()



