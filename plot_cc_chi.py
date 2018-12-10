import matplotlib.pyplot as plt
import sys
import math 
import numpy

datafile1 = sys.argv[1]

f = open(datafile1,'r')


cc =[]
chi_square_log = []

for ll in f.readlines():
  elements=ll.split()
  cc.append(float(elements[0]))
  chi_square_log.append(numpy.log(float(elements[0])))


plt.scatter(cc,chi_square_log)
plt.xlabel('cc') 
# plt.xlim(xmax=0.75,xmin=0.7)
plt.ylabel('chi_square_log') 


plt.savefig("./cc_chi_5484.png")
plt.show()