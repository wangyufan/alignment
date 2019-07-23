import numpy as np
from scipy.sparse import csr_matrix
import os, sys
import operator


if __name__=="__main__":
	args = sys.argv[1:]
	# arr = np.load("total_cc.npy")
	# row = [x[0] for x in arr]
	# col = [x[1] for x in arr]
	# data = [abs(1-abs(x[2])) for x in arr]
	# res1 = csr_matrix((data, (col, row)), shape=(5298, 5298))
	# res2 = csr_matrix((data, (row, col)), shape=(5298, 5298))
	# (res1.A +res2.A).tofile("yufan_distance.matrix") 
	
	# np.save("distance2.matrix", res.A)
	# np.set_printoptions(threshold=np.inf)
	# mat = np.load("distance2.matrix.npy")
	# mat2 = np.load("distance.matrix")
	# print(mat+mat2)

	# new_arr = [(x[0], x[1], abs(1-abs(x[2])))for x in arr]
	# cc_arr = [(x[0], x[1], abs(x[2]))for x in arr]
	# sort_cc_arr = sorted(cc_arr, key=operator.itemgetter(2), reverse=True)
	# np.savetxt("sort_cc_matrix.txt", sort_cc_arr, fmt="%d,%d,%f")
	# sort_arr = sorted(new_arr, key=operator.itemgetter(2))
	# np.savetxt("sort_distance.txt", sort_arr, fmt="%d,%d,%f")
	

	with open("yufan_distance.matrix", 'rb') as matrix:
		data = np.fromfile(matrix, count=5298*5298)
		array = np.reshape(data, [5298, 5298], order='F')
	print("Values bigger than 1 =", array[array > 1])
	print("Their indices are :", np.argwhere(array > 1))
	print("The matrix is : \n", array)
	# a = np.load("distance.matrix")
	# a2 = np.reshape(a, [5298, 5298], order='F')
	# print(a2)
	# print("Their indices are ", np.argwhere(a2 > 1))
