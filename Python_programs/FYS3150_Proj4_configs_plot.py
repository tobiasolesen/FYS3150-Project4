# -*- coding: utf-8 -*-

from __future__ import division
from matplotlib.pyplot import *
from numpy import *

config_list = []

count = 0

config_file = open("configs.txt", "r")

#Leser fra config_file:
for line in config_file:
    numbers = line.split()
    config_list.append(numbers[0])
    count +=1

config_file.close()


#config_list = array(config_list)

mcs = count
mcs_list = linspace(0, mcs, mcs)

print len(mcs_list)
print len(config_list)

plot(mcs_list[:1000], config_list[:1000])
title("Number of accepted configurations vs number of Monte Carlo cycles with an random starting configuration (T=2.4)")
xlabel("Number of Monte Carlo cycles")
ylabel("Number of accepted configurations")
show()