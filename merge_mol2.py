#coding:utf-8
import sys
reload(sys)
import os
import os.path
import time

time1=time.time()
def merge_file(filestr):
    targetfile = os.path.join(os.path.split(sys.path[0])[0],"mol2_merge.mol2") 
    print '=======================',targetfile
    tempf = open(targetfile,'w')
    filelist = filestr[0].split(",")
    print "filelist",filestr
    for i in range(len(filelist)):
    	f = open(filelist[i],'r')
    	tempf.write(f.read()+"\n")

	print u'merge    ' + str(i) + '  mol2 file finished'
    tempf.close()
    return 


if __name__ == '__main__':
    args=sys.argv[1:]
    merge_file(args)
    time2 = time.time()
    print u'time used' + str(time2 - time1) + 's'
 
