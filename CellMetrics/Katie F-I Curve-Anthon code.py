# -*- coding: utf-8 -*-
#Created on June 2, 2017
#Modified on June 2, 2017
#Author: Anton Lunyov
#Rewriting the code by Amir on Dec 2017

from brian2 import *
from pylab import *

#No compiling to speed up process
prefs.codegen.target = "numpy"
#Don't include previous Brian objects
start_scope()

# Strongly adapting params (Izhikevitch)
Cm = 90.0 * pF # Membrane capacitance
vR = -60.6 * mV # Resting membrane potential
vPeak = 2.5 * mV # Spike cutoff
vT = -43.1 * mV # Instantaneous membrane potential
Ishift = 0.0 * pA # Rheobase current shift
Iother = 0.0 * pA # Defined at zero since input is not noisy

kLow = 1.7 * nS/mV
kHigh = 14 * nS/mV

a = 0.1 /ms
d = 0.1 * pA
b = -0.1 * nS
c = -67 * mV

#Single neuron simulation - 60 individual neurons with varying current steps
N = 150
#Arrays for storing the initial and final ISI
freq_initial = zeros(N)
freq_final = zeros(N)

#The values of the default current input - max and min. Down-up-down
stepDown = 0.0 * pA
stepUp = 2.0 * pA

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

defaultclock.dt = 0.01*ms #Integration timecuts

#Set duration and run simulation
duration = 1.5 * second

run(duration)
count = 0 #counting the number of input currents for which the frequency is zero to find the rheobase current

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
    count += (freq_initial[h] < 1) * 1 #setting upper limit for rheobase current definition 

#Do a linear fit of the first 5 data points above 1 Hz
x = arange(0,stepUp*N/pA,stepUp/pA)
fitInitial = polyfit(x[count:count+5],freq_initial[count:count+5],1)
#fitFinal = polyfit(x,freq_final,1)
fitInitPlot = (fitInitial[0] * array(x)) + fitInitial[1]
#fitFinalPlot = (fitFinal[0] * array(x)) + fitFinal[1]

print('Cm = %s, ' %Cm, 'b= %s, ' %b, 'kLow= %s' %kLow ) #printing all the varying data
print('Irheo by definition = %s' %x[count]) #printing rheobase current using row data
print('Irheo by fitting line = %s' %((((-1)*fitInitial[1])/fitInitial[0]))) #printing rheobase current using graph

#Plot relevant data
figure(1) 
plot(x,freq_initial,'b', label='Model_initial freq')
plot(x,freq_final,'og', label='Model_final freq')
#plot(x,fitInitPlot,'g', label='Fitted line to first 5 data above 1 Hz')
xlim(xmin=0)
ylim(ymin=0)
legend()
xlabel("Input Current (pA)")
ylabel("Frequency (Hz)")
show()
