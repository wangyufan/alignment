#-*-coding:utf-8-*-

import os
from libtbx import easy_pickle
from scitbx.array_family import flex
import numpy as np
import matrix_pearson
import datetime

begin = datetime.datetime.now()
data = easy_pickle.load("//home/dongxq/align_code/dude-actives.nl")
length = len(data)

count = 0
diff_list = [ 0 for i in range(int(length*(length-1)/2))]

for i in range(length-1):
	for j in range(i+1,length):
		diff_listi = matrix_pearson.pearson_cc(data[i],data[j])
		diff_list[count] = diff_listi
		print diff_listi
		count = count + 1

# diff_array = np.array(diff_list)
#np.save('diffmol2-090.npy',diff_list)
end = datetime.datetime.now()
#print data
print str(end-begin)
print count
# print (data[0] - data[1]).norm()
print length