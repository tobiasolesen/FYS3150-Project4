# -*- coding: utf-8 -*-

from __future__ import division
from matplotlib.pyplot import *
#import matplotlib.pyplot as mplt
from numpy import *

E_file = open("output.txt", "r")

E_list = []
temp_list = []

#List of discrete E-values after equilibrium is reached (not divided by n_spins^2)
#For T=1:
#discrete_E = linspace(-720.,-800., 21) #21 punkter, gir steps of dE=4.Tar ikke med de forste energiene for starter ikke a beregne P(E) before steady state is reached

#For T = 2:
discrete_E = linspace(-540., -840., 76)

#For T = 2.4:
#discrete_E = linspace(-340., -640., 76) #Gir DE = 4

numb_discrete_E = len(discrete_E)  #Number of discrete E-values
E_count_vec = zeros(numb_discrete_E) #vektor som vil inneholde antallet av de forskjellige E-verdiene
prob_vec = zeros(numb_discrete_E) #vil inneholde sannsynlighetene for de forskjellige E-verdiene.


count = 0
n_spins = 20

#E_file.readline()

#Leser fra E_file:
for line in E_file:
    numbers = line.split()
    temp_list.append(float(numbers[0]))
    E_list.append(float(numbers[1]))
    count+=1
E_file.close()


#Counts number of times each of the energies occur in order to calculate the probability for a given energy
for i in range(len(E_list)):
    for j in range(numb_discrete_E):    
        if E_list[i] == discrete_E[j]: #For T = 2.4
            E_count_vec[j] += 1
        #if E_list[i] == discrete_E[j]: #For T = 1
            #E_count_vec[j] += 1
    
prob_vec = E_count_vec / sum(E_count_vec) #Vector with all probabilities
print "Sannsynlighetsvektor:", prob_vec
print len(prob_vec)
print len(discrete_E)
prob_list = list(prob_vec)
print type(prob_list)
 
#Sannsynlighet for E = -800 (som er den mest sannsynlige siden vi konvergerer mot den):
print "Sannsynlighet for E=", discrete_E[-1], ":", prob_vec[-1]
#plot(discrete_E, prob_vec, 'o') #For T=1
plot(discrete_E, prob_vec,'o') #For T=2.4
title("Probability of discrete energy-values for 20x20 spins and T=2")
xlabel("Discrete energies")
ylabel("Probability")

    
#Ser fra plot at energien vil konvergere mot -2, og deretter fluctuate rundt denne verdien.
#Med flere mc cycles kovergerer energien raskere?    
    

#mplt.show()
