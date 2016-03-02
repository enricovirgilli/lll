# -*- coding: iso-8859-15 -*-
import os, time, thread, numpy
import sys, getopt, os
from pylab import array, load, save, zeros, plot, savefig, close

EA_b4_bin = load(os.path.join("EA_lens.dat"))

len_EA = len(EA_b4_bin)

min_energy = 80

j = 0

EA_1 = []
Erange_1 = []
for i in range(1,min_energy):
    EA_1.append(0)

for i in range(min_energy,len_EA+min_energy):
    j = i - min_energy
    EA_1.append(EA_b4_bin[j, 1])
    Erange_1.append(EA_b4_bin[j,0])

for i in range(len_EA+min_energy,1+len_EA+min_energy+(len_EA+min_energy)/2):
    EA_1.append(0)

j_min = 0
j_max = 0
EA_binned = []

for i in range(min_energy,len_EA+min_energy):
    j_min = int(i - i/2)
    j_max = int(i + i/2)
    EA_binned.append(sum(EA_1[j_min:j_max+1]))


E1 = sum(EA_1[100:117])
E2 = sum(EA_1[118:137])
E3 = sum(EA_1[138:160])
E4 = sum(EA_1[161:187])
E5 = sum(EA_1[188:219])
E6 = sum(EA_1[220:256])
E7 = sum(EA_1[257:300])


print 108.5, 18, E1 
print 127.5, 20, E2 
print 149, 23, E3 
print 174, 27, E4 
print 203.5, 32, E5 
print 238, 37, E6 
print 278.5, 44, E7



a = len(EA_binned)
b = 0
EA_binned_petal = []
for i in range(0,a):
    b = EA_binned[i]/20.
    EA_binned_petal.append(b)

#for i in range(1,len_EA):
    #print Erange_1[i], EA_binned_petal[i] #for one petal.
    #print Erange_1[i], EA_binned[i] #for total lens.


    





