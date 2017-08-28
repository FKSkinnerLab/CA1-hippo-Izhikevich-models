#Created on June 26, 2017
#Modified on Aug 1, 2017
#Author: Anton Lunyov
#Pulls out specified data from 4D tensors generated for analyzing the effects of a,b,d,kLow on SFA, PIR and rheobase current

from pylab import *

res = 10

# a value iterations
aRes = 2.4E-4 # 1.0E-4
aVals = zeros(res)
# klow iterations
kRes = 0.02
kVals = zeros(res)
# b iters
bRes = 0.6
bVals = zeros(res)
# d iters
dRes = 2
dVals = zeros(res)

for s in range(res):
    aVals[s] = aRes * (s**2)
    bVals[s] = bRes * (s**2)
    kVals[s] = kRes * (s**2)
    dVals[s] = dRes * (s**2)


#Load all relevant data
transitionCurrents = load('TC_square_tensor.npy')
adaptation_rat = load('SFA_ratio_square_tensor.npy')
adaptation = load('SFA_square_tensor.npy')
rheo = load('rheo_square_tensor.npy')

#print(unique(adaptation))
#print(unique(rheo))
print(unique(transitionCurrents, return_counts = True))



#Redo properly given square tensor
###Parameters to pull out a single value from
##a = 0.0012
##b = 3.0
##d = 10.0
##kLow = 0.10
##
###Calculate the indeces in 4D tensor
##aInd = int((a/(aVals[res - 1] - aVals[0]))*len(aVals))
##bInd = int((b/(bVals[res - 1] - bVals[0]))*len(bVals))
##dInd = int((d/(dVals[res - 1] - dVals[0]))*len(dVals))
##kInd = int((kLow/(kVals[res - 1] - kVals[0]))*len(kVals))
##
##print(aInd)
##
##print(transitionCurrents[aInd,bInd,dInd,kInd])
##print(adaptation[aInd,bInd,dInd,kInd])
##print(adaptation_rat[aInd,bInd,dInd,kInd])
##print(rheo[aInd,bInd,dInd,kInd])


#std_SFA = std(adaptation)
#std_TC = std(transitionCurrents)
#std_rheo = std(rheo)

SFA_target = 0.46 #0.46 norm #5.86 rat
rheo_target = 4.0
TC_target = -6.5

SFA_wind = 0.45
rheo_wind = 3.0
tc_wind = 3.0

SFA_up = SFA_target + SFA_wind
SFA_low = SFA_target - SFA_wind

rheo_up = rheo_target + rheo_wind
rheo_low = rheo_target - rheo_wind

TC_up = TC_target + tc_wind
TC_low = TC_target - tc_wind

inRange = (transitionCurrents < TC_up)*(transitionCurrents > TC_low)
inRange = inRange * (adaptation < SFA_up)*(adaptation > SFA_low)
inRange = inRange * (rheo < rheo_up)*(rheo > rheo_low)

#Get the parameter values within
aValsWithin = aVals[nonzero(inRange)[0]]
bValsWithin = bVals[nonzero(inRange)[1]]
dValsWithin = dVals[nonzero(inRange)[2]]
kValsWithin = kVals[nonzero(inRange)[3]]

suptitle(str(float(len(aValsWithin))*100/10000) + "% of the models were explored out of all")
subplot(221)
hist(aValsWithin, bins = aVals)
xlabel('a')
ylabel('Frequency of parameter')

subplot(222)
hist(bValsWithin, bins = bVals)
xlabel('b')
ylabel('Frequency of parameter')

subplot(223)
hist(dValsWithin, bins = dVals)
xlabel('d')
ylabel('Frequency of parameter')

subplot(224)
hist(kValsWithin, bins = kVals)
xlabel('k')
ylabel('Frequency of parameter')

show()
