# -*- coding:utf-8 -*-

import numpy as np
from libtbx import easy_pickle

# from numpy import *

def pearson_cc(x,y):
	xlength = len(x)
	ylength = len(y)
	if xlength != ylength:
		print 'sample length is not equal!'
	# print xmean

	x = np.matrix(x)
	y = np.matrix(y)

	Sum_xy = np.sum( np.multiply((x-np.mean(x)),(y-np.mean(y))) )
	covariance = Sum_xy/(xlength-1)
	Sx = np.sum(np.multiply((x-np.mean(x)),(x-np.mean(x))))
	# print Sx
	stddx = pow(Sx/(xlength-1), 0.5)
	Sy = np.sum(np.multiply((y-np.mean(y)),(y-np.mean(y))))
	# print Sy
	stddy = pow(Sy/(ylength-1), 0.5)
	corr_xy = covariance/(stddx*stddy)

	# print corr_xy
	return corr_xy
if __name__ == '__main__':
	# x = [-3,6,0,3,-6]
	# y = [1,-2,0,-1,2]
	x = [3,6,0,3,6]
	y = [1,2,0,1,2]
	cc = pearson_cc(x,y)
	print cc

	# for i in range(xlength):
	# 	# print x[i]
	# 	Sum_xy += (x[i]-xmean)*(y[i]-ymean)
	# covariance = Sum_xy/(xlength-1)

	# for i in range(xlength):
	# 	Sx += (x[i]-xmean)*(x[i]-xmean)
	# stddx = pow(Sx/(xlength-1), 0.5)

	# for i in range(ylength):
	# 	Sy += (y[i]-xmean)*(y[i]-ymean)
	# stddy = pow(Sy/(ylength-1), 0.5)

	# corr_xy = covariance/(stddx*stddy)

	# return corr_xy