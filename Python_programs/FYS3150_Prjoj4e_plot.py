# -*- coding: utf-8 -*-

from __future__ import division
from matplotlib.pyplot import *
from numpy import *

#Alle verdier er per spin

#T_start = 2.0
#T_end = 2.2
#dT = 0.05

E_list = []
E2_list = []
M_list = []
M2_list = []
T_list = []
susceptibility_list = []
heatCapacity_list = []

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

file.close()

subplot(4, 1, 1)
title("Energy vs temperature (L = 100) for T=[2.1,2.4], dT = 0.05, mcs=1000000")
xlabel("Temperature")
ylabel("Energy per spin")
plot(T_list, E_list)

subplot(4, 1, 2)
title("Magnetization vs temperature")
plot(T_list, M_list)

subplot(4, 1, 3)
title("Heat capacity vs temperature")
plot(T_list, heatCapacity_list)

subplot(4, 1, 4)
title("Susceptibility vs temperature")
plot(T_list, susceptibility_list)

show()

