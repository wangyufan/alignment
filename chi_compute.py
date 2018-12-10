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

def get_mean_sigma( nlm_array ):
  coef = nlm_array.coefs()
  mean = abs( coef[0] )
  var = flex.sum( flex.norm(coef) )
  sigma = smath.sqrt( var-mean*mean )
  return mean, sigma

def computeChi(process_n):
#compute distance 
  # print "\nRun task chi-%s" %(os.getpid()), time.ctime()
  nlRes = []
  cavity_nl_array = cavity_model.nl_array
  cavity_nl_coefs = cavity_model.nl_array.coefs()
  start_line = process_n*perSizeChi
  if ((process_n+1)*perSizeChi < totalSize):  
    # range(x,y) => [x,y-1]
    end_line = (process_n+1)*perSizeChi     
  else:
    end_line = totalSize
  # print "chi_end_line:",end_line
  for indx in range(start_line, end_line):
  #compute Chi-sequare distance
    mf_coef = numpy.true_divide(nl_coefs[indx], cavity_nl_coefs)
    dist = numpy.sum(numpy.square(cavity_nl_coefs - mf_coef*nl_coefs[indx]))
  #compute Mahalanobis distance
    # dist = mol2.Mahalanobis(cavity_nl_coefs,nl_coefs[indx])
    nlRes.append((indx, dist, codes[indx]))
  return nlRes

def run(args):
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
  nlm_coefs=easy_pickle.load(prefix+".nlm")
  nl_coefs=easy_pickle.load(prefix+".nl")
  rmaxs=easy_pickle.load(prefix+".rmax")

  cavity = args.cavitypath
  cavitytype = args.cavitytype
  nmax = args.nmax
  rmax = args.rmax
  processnum = args.processnum
  cavity_model=model_interface.build_model( cavity, cavitytype, nmax, rmax )
  totalSize = len(nl_coefs)
  perSizeChi = int(xmath.ceil( totalSize / processnum))
  NlRes = []

#chi process pool
  tnl1 = time.time()
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
