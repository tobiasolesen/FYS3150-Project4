# -*- coding: utf-8 -*-

from __future__ import division
from matplotlib.pyplot import *
from numpy import *

E_file = open("output.txt", "r")
config_file = open("configs.txt", "r")

#4 counters for the different E-values because there are 4 differents E-values
count_0 = 0 #Teller E[0]
count_1 = 0 #Teller E[1]
count_2 = 0 #Teller E[2]
count_3 = 0 #Teller E[3]

E_list = []
#List of discrete E-values after equilibrium is reached (not divided by n_spins^2)
discrete_E = linspace(-720.,-800, 21) #Gir steps of dE=4.Tar ikke med de forste energiene for starter ikke a beregne P(E) before steady state is reached
numb_discrete_E = len(discrete_E)  #Number of discrete E-values
E_count_vec = zeros(numb_discrete_E) #vektor som vil inneholde antallet av de forskjellige E-verdiene
prob_vec = zeros(numb_discrete_E) #vil inneholde sannsynlighetene for de forskjellige E-verdiene.

M_list = []
config_list = []
count = 0
n_spins = 20

#Leser fra E_file:
for line in E_file:
    numbers = line.split()
    E_list.append(float(numbers[0]))
    M_list.append(float(numbers[1]))
    count+=1
E_file.close()

#Counts number of times each of the energies occur
for i in range(len(E_list)):
    for j in range(numb_discrete_E):    
        if E_list[i] == discrete_E[j]:
            E_count_vec[j] += 1
    
prob_vec = E_count_vec / sum(E_count_vec)
print "Sannsynlighetsvektor:", prob_vec

#print E_count_vec  / sum(E_count_vec) 
    
#total_numb_of_Evalues = count_0 + count_1 + count_2 + count_3 #number of E-values involved in our calculation

#Sannsynlighet for E = -800 (som er den mest sannsynlige siden vi konvergerer mot den):
print "Sannsynlighet for E=", discrete_E[-1], ":", prob_vec[-1]
    
#Leser fra config_file:
for line in config_file:
    numbers = line.split()
    config_list.append(numbers[0])
config_file.close()
   
mcs = count
mcs_list = linspace(0, mcs, mcs)
     
E = array(E_list)    
M = array(M_list) 
    
#Ser fra plot at energien vil konvergere mot -2, og deretter fluctuate rundt denne verdien.
#Med flere mc cycles kovergerer energien raskere?    
    
#plot(mcs_list, E/(n_spins**2)
plot(mcs_list, config_list)
xlabel("Number of monte carlo cycles")
ylabel("E (mean energy per spin)")
title("mean energy per spin vs number of mc-cycles (random start)")
show()