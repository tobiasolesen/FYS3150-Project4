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

M_list = []
E_list = []
#List of discrete E-values after equilibrium is reached (not divided by n_spins^2)
discrete_E = linspace(-720.,-800, 21) #Gir steps of dE=4.Tar ikke med de forste energiene for starter ikke a beregne P(E) before steady state is reached
numb_discrete_E = len(discrete_E)  #Number of discrete E-values
E_count_vec = zeros(numb_discrete_E) #vektor som vil inneholde antallet av de forskjellige E-verdiene
prob_vec = zeros(numb_discrete_E) #vil inneholde sannsynlighetene for de forskjellige E-verdiene.

T_values = [2.0, 2.05, 2.1, 2.15, 2.2] #dT = 0.05

temp_list = []
config_list = []
count = 0
n_spins = 40

#E_file.readline()

#For 

#Leser fra E_file:
for line in E_file:
    numbers = line.split()
    temp_list.append(float(numbers[0]))
    E_list.append(float(numbers[1]))
    #M_list.append(float(numbers[2])) #Gir feil!
    count+=1

E_file.close()

mcs = count/len(T_values) #Oppg e)
print "mcs:", mcs
#T_end = T_values[-1]
T_end = 2.5
T = linspace(0, T_end, mcs)


'''
#Counts number of times each of the energies occur in order to calculate the probability for a given energy
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
'''
mcs_list = linspace(0, mcs, mcs)
     
E = array(E_list)/(n_spins**2)
M = array(M_list) 
    
print len(E[2*mcs:3*mcs])
print len(T)
print len(E)
    
legend_skrift = "T ="
    
for i in range(len(T_values)):
    print i
    plot( T, E[(i)*mcs:((i+1)*mcs)], label= (legend_skrift, T_values[i]) ) #Plot energy per spin
    
#for T in T_values:
#    plot(T, E[:mcs]/(n_spins**2))    
    
#Ser fra plot at energien vil konvergere mot -2, og deretter fluctuate rundt denne verdien.
#Med flere mc cycles kovergerer energien raskere?    
    
#plot(mcs_list, E/(n_spins**2)
#plot(mcs_list, M)
#plot(mcs_list, config_list)
#plot(mcs_list, E/(n_spins**2))
legend()
xlabel("Temperature")
#xlabel("Number of monte carlo cycles")
ylabel("E (mean energy per spin)")
#title("mean energy per spin vs number of mc-cycles (random start)")
title("Mean energy per spin vs temperature (random start)")
show()
