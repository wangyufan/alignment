from libtbx import easy_pickle
import os, sys
import model_interface
import mol22zernike
from scitbx import math, matrix
from iotbx import pdb
import mol2
import iotbx.phil
from scitbx.math import zernike_align_fft as fft_align
import time
import operator
import numpy
import math as xmath
import multiprocessing
from multiprocessing import Pool
# import pathos.multiprocessing 
# from pathos.multiprocessing import ProcessingPool as pPool
import argparse


banner = "-------------------Align molecular models in various format (PDB/mol2/NLM) ---------------------"


def computeChi(process_n):
#compute distance 
  # print "\nRun task chi-%s" %(os.getpid()), time.ctime()

  # cavity_model, processnum, process_n
  nlRes = []
  cavity_rmax = cavity_model.rmax
  totalSize = len(rmaxs)
  perSizeChi = int(xmath.ceil( totalSize / processnum))
  start_line = process_n*perSizeChi
  if ((process_n+1)*perSizeChi < totalSize):  
    # range(x,y) => [x,y-1]
    end_line = (process_n+1)*perSizeChi     
  else:
    end_line = totalSize
  # print "chi_end_line:",end_line
  rmax_dict = {}
  for i in range(totalSize):
      rmax_dict[i] = rmaxs[i]
  # print "----",end_line
  for indx in range(start_line, end_line):
    if indx == 25121:
      print "rmax in database:",rmaxs[indx]
      print "rmax of cavity:",cavity_rmax*0.9
      print "dist:",round(abs(rmaxs[indx] - cavity_rmax*0.9), 6)
    dist = round(abs(rmaxs[indx] - cavity_rmax*0.9), 6)
    nlRes.append((indx, dist, codes[indx]))
  sortedNlRes=sorted(nlRes, key=operator.itemgetter(1), reverse=False)
  return sortedNlRes

def run(args):
  global processnum, cavity_model, rmaxs, codes
  t1 = time.time()
  parser = argparse.ArgumentParser()
  parser.add_argument("-output", help="chi result path", type=str)
  parser.add_argument("-dbpath", help="db path", type=str)
  parser.add_argument("-processnum", help="processnum number", type=int)
  parser.add_argument("-cavitypath", help="cavity path", type=str)
  parser.add_argument("-cavitytype", help="cavity type", type=str)
  parser.add_argument("--nmax", help="nmax", default=20, type=int)
  parser.add_argument("--rmax", help="rmax", type=float)
  args = parser.parse_args()


  # prefix="/home/dongxq/zalign/build/myDB"
  # prefix="/home/dongxq/align_code/myDB/dude2-db"
  targetfile = args.output
  prefix = args.dbpath 
  codes=easy_pickle.load(prefix+".codes")
  rmaxs=easy_pickle.load(prefix+".rmax")
  cavity = args.cavitypath
  cavitytype = args.cavitytype
  nmax = args.nmax
  rmax = args.rmax
  processnum = args.processnum
  cavity_model=model_interface.build_model( cavity, cavitytype, nmax, rmax )
  NlRes = []

#chi process pool
  tnl1 = time.time()
  # computeChi(rmaxs, cavity_model, processnum)
  pool_chi = multiprocessing.Pool(processes = processnum)
  chi_result = pool_chi.map(computeChi, range(processnum))
  for n in range(0,len(chi_result)):
    NlRes += chi_result[n]
  tnl2 = time.time()    

#output
  with open(targetfile,"w") as f:
    for arr in NlRes:
      f.write(str(arr) + "\n")
    t2 = time.time()
    # f.write( "rotation invariant computing time used: "+str(tnl2-tnl1)+"\n")
    # f.write( "total time used: : "+str(t2-t1))
if __name__=="__main__":
  args = sys.argv[1:]
  run(args)
