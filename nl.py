from libtbx import easy_pickle
import os, sys
import model_interface
import mol22zernike
import mol2
import iotbx.phil
from scitbx import math, matrix
from scitbx.math import zernike_align_fft as fft_align
import time
import operator
import numpy
import math as xmath
import multiprocessing
from multiprocessing import Pool
import argparse

banner = "-------------------Align molecular models in various format (PDB/mol2/NLM) ---------------------"

def computeChi(process_n):
#compute distance 
  # print "\nRun task chi-%s" %(os.getpid()), time.ctime()
  nlRes = []
  mov_nl_coefs = mov_model.nl_array.coefs()
  start_line = process_n*perSizeChi
  if ((process_n+1)*perSizeChi < totalSize):  
    # range(x,y) => [x,y-1]
    end_line = (process_n+1)*perSizeChi     
  else:
    end_line = totalSize
  # print "chi_end_line:",end_line
  for indx in range(start_line, end_line):
  #compute Chi-sequare distance
    mf_coef = numpy.true_divide(nl_coefs[indx], mov_nl_coefs)
    dist = numpy.sum(numpy.square(mov_nl_coefs - mf_coef*nl_coefs[indx]))
  #compute Mahalanobis distance
    # dist = mol2.Mahalanobis(mov_nl_coefs,nl_coefs[indx])
    nlRes.append((indx, dist, codes[indx]))
  return nlRes



def computeCc(process_n):
#compute nlm_cc
  print "\nRun task cc pid%s" %(os.getpid()), time.ctime() 
  mov_nlm_array = mov_model.nlm_array
  fix_nlm_array = math.nlm_array(nmax)
  nlm = fix_nlm_array.nlm()
  nlm_total = fix_nlm_array.nlm().size()
  nlmRes = []
  
  start_line = process_n*perSizeCc
  if (process_n+1)*perSizeCc > totalSize:
    end_line = totalSize
  else:
    end_line = (process_n+1)*perSizeCc
  print "cc_end_line:",end_line
  for i in range(start_line, end_line):
    indx = sortedNlRes[i][0]
    fix = nlm_coefs[indx][0:nlm_total]
    fix_nlm_array.load_coefs(nlm, fix)
    align_obj = fft_align.align( fix_nlm_array, mov_nlm_array, nmax=nmax, refine=True )
    cc = align_obj.get_cc()   
    nlmRes.append((indx, codes[indx], cc))
  # print "thread %d finished at"%(process_n), time.ctime()
  return nlmRes

if __name__=="__main__":
  t1 = time.time()

  parser = argparse.ArgumentParser()
  parser.add_argument("-output", help="cc result path", type=str)
  # parser.add_argument("-num_fnl_round", help="nl", type=int)
  # parser.add_argument("-num_nn_round", help="nn number", type=int)
  parser.add_argument("-processnum", help="processnum number", type=int)
  parser.add_argument("-dbpath", help="db path", type=str)
  parser.add_argument("-target_shape", help="cavity path", type=str)
  parser.add_argument("-target_type", help="cavity type", type=str)
  parser.add_argument("--nmax", help="nmax", default=20, type=int)
  parser.add_argument("--rmax", help="rmax", type=float)
  args = parser.parse_args()
  
  # prefix="/home/dongxq/zalign/build/myDB"
  # prefix="/home/wyf/project/align_code/myDB/dude-actives"
  prefix=args.dbpath
  db_name = os.path.basename(prefix)
  codes=easy_pickle.load(prefix+".codes")
  nlm_coefs=easy_pickle.load(prefix+".nlm")
  nl_coefs=easy_pickle.load(prefix+".nl")
  rmaxs=easy_pickle.load(prefix+".rmax")
  totalSize = len(nl_coefs)
  nmax = args.nmax
  rmax = args.rmax
  processnum = args.processnum
  perSizeChi = int(xmath.ceil( totalSize / processnum))
  perSizeCc = int(xmath.ceil(totalSize / processnum))
  mov_model=model_interface.build_model( args.target_shape, args.target_type, nmax, rmax )
  sortedNlRes = []
  sortedNlmRes = []


#chi process pool
  tnl1 = time.time()
  pool_chi = multiprocessing.Pool(processes = processnum)
  chi_result = pool_chi.map(computeChi, range(processnum))
  for n in range(0,len(chi_result)):
    sortedNlRes += chi_result[n]
  # sortedNlRes=sorted(sortedNlRes, key=operator.itemgetter(1), reverse=False)[:args.num_fnl_round]
  tnl2 = time.time()    

#cc process pool
  tnlm1 = time.time()
  pool_cc = multiprocessing.Pool(processes=processnum)
  cc_result = pool_cc.map(computeCc, range(processnum))
  for n in range(0,len(cc_result)):
    sortedNlmRes += cc_result[n]
  # sortedNlmRes = sorted(sortedNlmRes, key=operator.itemgetter(2), reverse=True)[:args.num_nn_round]
  tnlm2 = time.time()

#merge chi to cc arr
  tmerge1 = time.time()
  for i in range(totalSize):
    indx = sortedNlmRes[i][0]
    chi = list(filter(lambda j: j[0] == indx, sortedNlRes[0:]))[0][1]
    sortedNlmRes[i] += (chi,)
  tmerge2 = time.time()

#output
  with open(args.output,"w") as f:
    rank = 0
    for arr in sortedNlmRes:
      f.write(str(arr) + "\n")
    t2 = time.time()
  print "rotation invariant computing time used: ", db_name, str(tnl2-tnl1)+"\n"
  print "alignment computing time used: ", db_name, str(tnlm2-tnlm1)+"\n"
  print "merge time used: " ,db_name, str(tmerge2-tmerge1)+"\n"
  print "total time used: : ",db_name,  str(t2-t1)