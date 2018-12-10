import matplotlib.pyplot as plt
import sys
import math 
import numpy

datafile1 = sys.argv[1]
datafile2 = sys.argv[2]

f1 = open(datafile1,'r')
f2 = open(datafile2,'r')

# cc = 0
# avg_cc = 0
# avg_cc_arr =[]
# model = 0
# model_arr = []
# cc2 = 0
# avg_cc2 = 0
# avg_cc_arr2 =[]
# model2 = 0
# model_arr2 = []

def cc_model_plot(f, select):
	cc = 0
	avg_cc = 0
	avg_cc_arr =[]
	model = 0
	model_arr = []
	for ll in f.readlines():
	  model +=1
	  elements=ll.split()
	  cc +=float(elements[select])
	  avg_cc = numpy.true_divide(cc, model)
	  avg_cc_arr.append(avg_cc)
	  model_arr.append(model)
	  # if model == 200:
	  # 	break
	return avg_cc_arr, model_arr

avg_cc_arr1, model_arr1 = cc_model_plot(f1,0)
avg_cc_arr2, model_arr2 = cc_model_plot(f2,0)


# for ll in f1.readlines():
#   model +=1
#   elements=ll.split()
#   cc +=float(elements[0])
#   avg_cc = numpy.true_divide(cc, model)
#   avg_cc_arr.append(avg_cc)
#   model_arr.append(model)
#   if model == 200:
#   	break

# for ll in f2.readlines():
#   model2 +=1
#   elements=ll.split()
#   cc2 +=float(elements[0])
#   avg_cc2 = numpy.true_divide(cc2, model2)
#   avg_cc_arr2.append(avg_cc2)
#   model_arr2.append(model2)
#   if model2 == 200:
#   	break

# plt.scatter(cc,chi_square_log)
plt.xlabel('top models') 
# plt.xlim(xmax=0.75,xmin=0.7)
plt.ylabel('avg_cc') 

#point
plt.scatter(model_arr1, avg_cc_arr1, s=20)
plt.scatter(model_arr2, avg_cc_arr2, s=20)

#curve 
# plt.plot(model_arr1, avg_cc_arr1,label='pdb5000_cc')
# plt.plot(model_arr2, avg_cc_arr2,label='pdb200_chi-square')
# plt.legend() 
plt.savefig("./avg_cc_points_5000.png")
plt.show()