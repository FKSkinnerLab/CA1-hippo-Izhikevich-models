#Created on June 26, 2017
#Modified on June 26, 2017
#Author: Anton Lunyov
#Pulls out specified data from 4D tensors generated for analyzing the effects of a,b,d,kLow on SFA, PIR and rheobase current

from pylab import *


##obtained from bubble_rheo_tensor.py, but the same for all (for consistency)
# a value iterations
aRes = 2.4E-4 # 1.0E-4
aUp = 2.4E-3
aLow = 0.0E-3
aVals = arange(aLow,aUp,aRes)

# klow iterations
kRes = 0.02 # 0.01
kLowUp = 0.2
kLowLower = 0.0
kVals = arange(kLowLower,kLowUp,kRes)

# b iters
bRes = 0.6 #0.2
bUp = 6.0
bLow = 0.0
bVals = arange(bLow,bUp,bRes)

# d iters
dRes = 2 # 0.5
dUp = 20
dLow = 0
dVals = arange(dLow,dUp,dRes)


#Load all relevant data
transitionCurrents = load('TC_across_4.npy')
adaptation_rat = load('adaptation_across_4_ratio.npy')
adaptation = load('adaptation_across_4.npy')
rheo = load('rheo_across_4.npy')

#Parameters to pull out a single value from
a = 0.0012
b = 3.0
d = 10.0
kLow = 0.10

#Calculate the indeces in 4D tensor
aInd = int((a/(aUp - aLow))*len(aVals)) - 1
bInd = int((b/(bUp - bLow))*len(bVals)) - 1
dInd = int((d/(dUp - dLow))*len(dVals)) - 1
kInd = int((kLow/(kLowUp - kLowLower))*len(kVals)) - 1

print(transitionCurrents[aInd,bInd,dInd,kInd])
print(adaptation[aInd,bInd,dInd,kInd])
print(adaptation_rat[aInd,bInd,dInd,kInd])
print(rheo[aInd,bInd,dInd,kInd])


std_SFA = std(adaptation)
std_TC = std(transitionCurrents)
std_rheo = std(rheo)

SFA_target = 0.46 #0.46 norm #5.96 rat
rheo_target = 3.5
TC_target = -6.0

SFA_up = SFA_target + 0.5*std_SFA
SFA_low = SFA_target - 0.5*std_SFA

rheo_up = rheo_target + 0.5*std_rheo
rheo_low = rheo_target - 0.5*std_rheo

TC_up = TC_target + 0.5*std_TC
TC_low = TC_target - 0.5*std_TC

inRange = (transitionCurrents < TC_up)*(transitionCurrents > TC_low)
inRange = inRange * (adaptation < SFA_up)*(adaptation > SFA_low)
inRange = inRange * (rheo < rheo_up)*(rheo > rheo_low)
#inRange = inRange*adaptation

#Get the parameter values within
aValsWithin = aVals[nonzero(inRange)[0]]
bValsWithin = bVals[nonzero(inRange)[1]]
dValsWithin = dVals[nonzero(inRange)[2]]
kValsWithin = kVals[nonzero(inRange)[3]]

suptitle(str(float(len(aValsWithin))*100/10000) + "% of the models were explored out of all")
subplot(221)
hist(aValsWithin)
xlabel('a')
ylabel('Frequency of parameter')

subplot(222)
hist(bValsWithin)
xlabel('b')
ylabel('Frequency of parameter')

subplot(223)
hist(dValsWithin)
xlabel('d')
ylabel('Frequency of parameter')

subplot(224)
hist(kValsWithin)
xlabel('k')
ylabel('Frequency of parameter')

show()
