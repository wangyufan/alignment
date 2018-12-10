import matplotlib.pyplot as plt
import sys
import math 

datafile = sys.argv[1]
plot_type = sys.argv[2]

f = open(datafile,'r')
print datafile
chi_square=[]
chi_square_log=[]
align_cc = []
nl_cc = []

for ll in f.readlines():
  elements=ll.split()
  align_cc.append(float(elements[0]))
  nl_cc.append(float(elements[1]))
  if float(elements[1]) > float(0.0):
  	print math.log(float(elements[1]))
  	chi_square_log.append(math.log(float(elements[1])))	
  else:
  	print float(elements[1])
  	chi_square_log.append(-15)

# plt.scatter(cc,chi_square_log)
plt.xlabel('align_cc') 
# plt.xlim(xmin=0.65,xmax=0.75)
if plot_type == '2':  
  plt.ylabel('nl_cc') 
  plt.scatter(align_cc,nl_cc)
else:
  plt.ylabel('log(Chi-square)') 
  plt.scatter(align_cc,chi_square_log)
filename = datafile.split('\n')[0].split('.')[0]
plt.savefig("./png"+str(filename)+".png")
plt.show()