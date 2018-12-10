#coding:utf-8
import sys
reload(sys)
import os
import os.path
import time
import argparse
import pickle
import operator



time1=time.time()
def file_sort(args):
	nlNum = args.nlnum
	inputnum = args.inputnum
	chiFileDir = args.input
	reverseFlag = args.reverse
	nlRes=[]
	#dblist=["dude1-db","dude2-db","dude-actives"]
	dblist=["anti3","antinon3"]
	for i in range(inputnum):
		chi_file = []
		chi_path = chiFileDir +"_"+ dblist[i]
		for line in open(chi_path, 'r').readlines():
			line = line.strip('\n')
			chi_file.append(eval(line) + (dblist[i],))	
		nlRes.extend(chi_file)
	if reverseFlag:
		sortKey = 2
	else:
	    sortKey = 1
	sortedNlRes=sorted(nlRes, key=operator.itemgetter(sortKey), reverse=reverseFlag)[:nlNum]
	# sortedNlRes = nlRes.sort(key=takeSecond, reverse=False)
	print args.reverse
	return sortedNlRes


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument("-output", help="sorted chi/cc path", type=str)
	parser.add_argument("-input", help="chi/cc file paths' root path", type=str)
	parser.add_argument("-inputnum", help="chi/cc file num", type=int)
	parser.add_argument("-nlnum", help="nlnum", type=int)
	parser.add_argument("--reverse", default=False, help=" a flag means to reverse the result, no matter T/F,the value is True", type=bool)
	# parser.add_argument('-input', action='append', dest='collection',
	#                 default=[],
	#                 help='Add chi file paths to a list',
	#                 )
	args = parser.parse_args()
	sortedNlRes = file_sort(args)
	targetfile = args.output
	with open(targetfile,"w") as f:
	    rank = 0
	    for arr in sortedNlRes:
	      rank +=1
	      arr=(rank,)+arr
	      f.write(str(arr) + "\n")
	time2 = time.time()
