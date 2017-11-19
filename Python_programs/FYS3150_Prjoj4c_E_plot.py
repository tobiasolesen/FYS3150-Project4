# -*- coding: utf-8 -*-

from __future__ import division
from matplotlib.pyplot import *
from numpy import *

E_list = []
E2_list = []
M_list = []
M2_list = []
T_list = []
susceptibility_list = []
heatCapacity_list = []

count = 0

file = open("output.txt", "r")

for line in file:
    numbers = line.split()
    T_list.append(float(numbers[0]))
    E_list.append(float(numbers[1]))
    E2_list.append(float(numbers[2]))
    M_list.append(float(numbers[3]))
    M2_list.append(float(numbers[4]))
    heatCapacity_list.append(float(numbers[5]))
    susceptibility_list.append(float(numbers[6]))
    count += 1

file.close()

mcs = count
mcs_list = linspace(0, mcs, mcs)

print len(mcs_list)
print len(E_list)

plot(mcs_list, E_list, label="E")
plot(mcs_list, M_list, label="M")

title("Energy and magnetization (per spin) vs number of monte carlo cycles, for T=2.4")
xlabel("Number of monte carlo cycles")
ylabel("Energy, magnetization")
legend()

show()