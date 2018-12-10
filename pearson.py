# -*- coding:utf-8 -*-

import datetime
from libtbx import easy_pickle
from scitbx.array_family import flex

def pearson_cc(x,y):
	xlength = len(x)
	ylength = len(y)
	# print x
	# print y
	if xlength != ylength:
		print 'sample length is not equal!'
	xmean = sum(x)/xlength
	ymean = sum(y)/ylength
	# print xmean

	Sum_xy = 0
	Sx = 0
	Sy = 0


	for i in range(xlength):
		# print x[i]
		Sum_xy += (x[i]-xmean)*(y[i]-ymean)
	covariance = Sum_xy/(xlength-1)

	for i in range(xlength):
		Sx += (x[i]-xmean)*(x[i]-xmean)
	stddx = pow(Sx/(xlength-1), 0.5)

	for i in range(ylength):
		Sy += (y[i]-xmean)*(y[i]-ymean)
	stddy = pow(Sy/(ylength-1), 0.5)

	corr_xy = covariance/(stddx*stddy)

	return corr_xy

# begin = datetime.datetime.now()
# data = easy_pickle.load("~/Desktop/mol2db/xiao_chwin.rmax")
# length = len(data)
# print dir(data)
# data_type = type(data)
# for i in range(length):
# 	print data[i]

# x = [-3,6,0,3,-6]
# y = [1,-2,0,-1,2]

# cc = pearson_cc(x,y)
# print cc
# print max(x)/3
# print pow(max(x),0.5)