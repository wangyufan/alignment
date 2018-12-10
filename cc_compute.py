from libtbx import easy_pickle
from stdlib import math as smath
from scitbx.array_family import flex
import os, sys
# from sastbx.zernike_model import model_interface
import model_interface
import mol22zernike

from scitbx import math, matrix
from iotbx import pdb
import mol2
import pearson
from sastbx.zernike_model import build_pymol_script

import iotbx.phil
from scitbx.math import zernike_align_fft as fft_align
from iotbx.xplor import map as xplor_map
from cctbx import uctbx
from scitbx.golden_section_search import gss
import time
import operator
import matplotlib
import matplotlib.pyplot as plt
import numpy
import math as xmath
import multiprocessing
from multiprocessing import Pool
import argparse


banner = "-------------------Align molecular models in various format (PDB/mol2/NLM) ---------------------"

def computeCc(process_n):
#compute nlm_cc
  print "\nRun task cc pid%s" %(os.getpid()), time.ctime() 
  cavity_nlm_array = cavity_model.nlm_array
  fix_nlm_array = math.nlm_array(nmax)
  nlm = fix_nlm_array.nlm()
  nlm_total = fix_nlm_array.nlm().size()
  nlmRes = []
  
  perSizeCc = int(xmath.ceil(nlNum / processnum))
  print "persize", perSizeCc
  start_line = process_n*perSizeCc
  if (process_n+1)*perSizeCc > nlNum:
    end_line = nlNum
  else:
    end_line = (process_n+1)*perSizeCc
  print "cc_end_line:",end_line
  for i in range(start_line, end_line):
    indx = sortedNlRes[i][1]
    fix = nlm_coefs[indx][0:nlm_total]
    fix_nlm_array.load_coefs(nlm, fix)
    align_obj = fft_align.align( fix_nlm_array, cavity_nlm_array, nmax=nmax, refine=True )
    cc = align_obj.get_cc()
    nlmRes.append((indx, codes[indx], cc))
  # print "thread %d finished at"%(process_n), time.ctime()
  return nlmRes

def run(args):
  t1 = time.time()
  # parser = argparse.ArgumentParser()
  # parser.add_argument("-output", help="cc result path", type=str)
  # parser.add_argument("-chipath", help="sorted chi path", type=str)
  # parser.add_argument("-node", help="node number", type=int)
  # parser.add_argument("-processnum", help="processnum number", type=int)
  # parser.add_argument("-dbpath", help="db path", type=str)
  # parser.add_argument("-cavitypath", help="cavity path", type=str)
  # parser.add_argument("-cavitytype", help="cavity type", type=str)
  # parser.add_argument("--nmax", help="nmax", default=20, type=int)
  # parser.add_argument("--rmax", help="rmax", type=float)
  # args = parser.parse_args()

  # prefix="/home/dongxq/align_code/myDB/dude1-db"
  prefix=args.dbpath
  db_name = os.path.basename(prefix)

  codes=easy_pickle.load(prefix+".codes")
  nlm_coefs=easy_pickle.load(prefix+".nlm")
  nl_coefs=easy_pickle.load(prefix+".nl")
  rmaxs=easy_pickle.load(prefix+".rmax")

  cavity = args.cavitypath
  cavitytype = args.cavitytype
  chipath = args.chipath
  nmax = args.nmax
  rmax = args.rmax
  processnum = args.processnum
  

  sortedNlRes = []
  countline = 0
  for line in open(chipath, 'r').readlines():
    line = line.strip('\n')
    line_tuple = eval(line)
    if line_tuple[-1] == db_name:
      countline += 1
      #can't del value in the tuple,so new one
      new_line = (countline,) + line_tuple[1:]
      sortedNlRes.append(new_line)       
  nlNum = len(sortedNlRes)
  print "------------",nlNum
  cavity_model=model_interface.build_model( cavity, cavitytype, nmax, rmax )
  sortedNlmRes = []

#cc process pool
  tnlm1 = time.time()
  pool_cc = multiprocessing.Pool(processes=processnum)
  cc_result = pool_cc.map(computeCc, range(processnum))

  for n in range(0,len(cc_result)):
    sortedNlmRes += cc_result[n]
  tnlm2 = time.time()

#merge chi to cc arr
  tmerge1 = time.time()
  for i in range(nlNum):
    indx = sortedNlmRes[i][0]
    chi = list(filter(lambda j: j[1] == indx, sortedNlRes[0:]))[0][2]
    sortedNlmRes[i] += (chi,)
  tmerge2 = time.time()

#output
  targetfile = args.output
  with open(targetfile,"w") as f:
    for arr in sortedNlmRes:
      f.write(str(arr) + "\n")
    t2 = time.time()

if __name__=="__main__":
  args = sys.argv[1:]
  run(args)
