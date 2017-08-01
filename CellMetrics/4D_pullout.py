#Created on June 26, 2017
#Modified on Aug 1, 2017
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

#Load all relevant data. Absolute values for adaptation can be obtained from adaptation_across_4.npy
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
aInd = int((a/(aUp - aLow))*len(aVals))
bInd = int((b/(bUp - bLow))*len(bVals))
dInd = int((d/(dUp - dLow))*len(dVals))
kInd = int((kLow/(kLowUp - kLowLower))*len(kVals))

#Output of relevant data
print("a: " + str(a))
print("b: " + str(b))
print("d: " + str(d))
print("kLow: " + str(kLow))

print("SFA: " + str(adaptation[aInd,bInd,dInd,kInd]))
print("PIR: " + str(transitionCurrents[aInd,bInd,dInd,kInd]))
print("Rheobase current: " + str(rheo[aInd,bInd,dInd,kInd]))


#Plots a specified 1D pullout, in this case k
subplot(311)
plot(dVals,adaptation[aInd,bInd,:,kInd],'go')
ylabel('SFA')

subplot(312)
plot(dVals,adaptation_rat[aInd,bInd,:,kInd],'go')
ylabel('SFA ratio')

subplot(313)
plot(dVals,rheo[aInd,bInd,:,kInd],'bo')
xlabel('d')
ylabel('Rheobase current')
show()
