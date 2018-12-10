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
from sastbx.interface import get_input
import matplotlib
import matplotlib.pyplot as plt
import numpy
import threading
import math as xmath




banner = "-------------------Align molecular models in various format (PDB/mol2/NLM) ---------------------"

def get_mean_sigma( nlm_array ):
  coef = nlm_array.coefs()
  mean = abs( coef[0] )
  var = flex.sum( flex.norm(coef) )
  sigma = smath.sqrt( var-mean*mean )
  return mean, sigma

def computeChi(thread_n, arr_chi):
#compute distance 
  print"============thread_chi============",thread_n
  nlRes = []
  mov_nl_array = mov_model.nl_array
  mov_nl_coefs = mov_model.nl_array.coefs()

  start_line = thread_n*perSizeChi
  if ((thread_n+1)*perSizeChi < totalSize):  
    # range(x,y) => [x,y-1]
    end_line = (thread_n+1)*perSizeChi + 1     
  else:
    end_line = totalSize

  for indx in range(start_line, end_line):
  #compute Chi-sequare distance
    mf_coef = numpy.true_divide(nl_coefs[indx], mov_nl_coefs)
    dist = numpy.sum(numpy.square(mov_nl_coefs - mf_coef*nl_coefs[indx]))
  #compute Mahalanobis distance
    # dist = mol2.Mahalanobis(mov_nl_coefs,nl_coefs[indx])
    nlRes.append((indx, dist, codes[indx]))
  arr_chi.extend(nlRes)


def computeCc(thread_n, arr_cc):
#compute nlm_cc
  print"----------thread_cc----------",thread_n
  mov_nlm_array = mov_model.nlm_array
  fix_nlm_array = math.nlm_array(nmax)
  nlm = fix_nlm_array.nlm()
  nlm_total = fix_nlm_array.nlm().size()
  nlmRes = []

  if ((thread_n+1)*perSizeCc < nlNum):  
    # range(x,y) => [x,y-1]
    end_line = (thread_n+1)*perSizeCc     
  else:
    end_line = nlNum

  for i in range(start_line, end_line):
    indx = sortedNlRes[i][0]
    fix = nlm_coefs[indx][0:nlm_total]
    fix_nlm_array.load_coefs(nlm, fix)
    align_obj = fft_align.align( fix_nlm_array, mov_nlm_array, nmax=nmax, refine=True )
    cc = align_obj.get_cc()
    nlmRes.append((indx, codes[indx], cc))
  arr_cc.extend(nlmRes)


if __name__=="__main__":
  t1 = time.time()
  base_path = os.path.split(sys.path[0])[0]
  global targetfile
  master_params = iotbx.phil.parse("""\
  align{
    nlnum = 10
    .type=int
    .help="nlNumnlNumnlNum"

    nlmnum = 10
    .type=int
    .help="nlmNum--"

    threadnum = 10
    .type=int
    .help="threadnum--"
     
    fix = None
    .type=path
    .help="pickle Cnlm coef of the fixed object"

    mov = None
    .type=path
    .help="pickle Cnlm coef of the fixed object"

    typef=*pdb mol2 nlm map
    .type=choice
    .help="fixed model type, default PDB"

    typem=*pdb mol2 nlm map 
    .type=choice
    .help="moving model type, default PDB"

    rmax = None
    .type=float
    .help="maxium radial distance to the C.o.M. (before the scaling)"

    nmax = 20
    .type=int
    .help="maximum order of zernike polynomial:fixed for the existing database"

    topn = 10
    .type=int
    .help="top N alignments will be further refined if required"

    refine = True
    .type=bool
    .help="Refine the initial alignments or not"

    write_map = False
    .type=bool
    .help="write xplor map to file"

  }
  output = "output"
  .type=path
  .help = "Output base. expect a .pr and a .qii file."
  """)

  args = sys.argv[1:]
  prefix="/home/dongxq/zalign/build/myDB"
  codes=easy_pickle.load(prefix+".codes")
  nlm_coefs=easy_pickle.load(prefix+".nlm")
  nl_coefs=easy_pickle.load(prefix+".nl")
  rmaxs=easy_pickle.load(prefix+".rmax")

  targetfile = os.path.join(os.path.split(sys.path[0])[0],"res")
  with open(targetfile,"w") as f:
    f.truncate()
  tempf = open(targetfile,'w')
  params = get_input(args, master_params, "aligndb", banner, help,tempf)
  tempf.close()

  mov = params.align.mov
  typem = params.align.typem
  nmax = params.align.nmax
  rmax = params.align.rmax
  topn = params.align.topn
  write_map = params.align.write_map
  nlNum = params.align.nlnum
  nlmNum = params.align.nlmnum
  threadNum = params.align.threadnum
  mov_model=model_interface.build_model( mov, typem, nmax, rmax )
  totalSize = len(nl_coefs)
  perSizeChi = int(xmath.ceil( totalSize / threadNum))
  perSizeCc = int(xmath.ceil(nlNum / threadNum))
  chi_arr = []
  cc_arr = []
  threads_chi = []
  threads_cc = []

#chi threads
  tnl1 = time.time()
  for i in range(int(threadNum)):
      threads_chi.append(threading.Thread(target=computeChi, args=(i, chi_arr,)))
  for t in threads_chi:
      t.setDaemon(True)
      t.start()
  for t in threads_chi:
      t.join() 

  sortedNlRes=sorted(chi_arr, key=operator.itemgetter(1), reverse=False)[:nlNum]
  tnl2 = time.time()

#cc threads
  tnlm1 = time.time()
  for j in range(int(threadNum)):
      threads_cc.append(threading.Thread(target=computeCc, args=(j, cc_arr,)))
  for h in threads_cc:
      h.setDaemon(True)
      h.start()
  for h in threads_cc:
      h.join()

  sortedNlmRes = sorted(cc_arr, key=operator.itemgetter(2), reverse=True)[:nlmNum]
  tnlm2 = time.time()
  print cc_arr
#merge chi to cc arr
  tmerge1 = time.time()
  for i in range(nlmNum):
    indx = sortedNlmRes[i][0]
    chi = list(filter(lambda j: j[0] == indx, sortedNlRes[0:]))[0][1]
    sortedNlmRes[i] += (chi,)
  tmerge2 = time.time()

#output
  with open(targetfile,"w") as f:
    f.write( "#############  SUMMARY of ALIGNMENT  #############\n")
    f.write( "rank indx     name               cc                 chi-square\n")
    rank = 0
    for arr in sortedNlmRes:
      rank +=1
      arr=(rank,)+arr
      f.write(str(arr) + "\n")
    t2 = time.time()
    f.write( "rotation invariant computing time used: "+str(tnl2-tnl1)+"\n")
    f.write( "alignment computing time used: "+str(tnlm2-tnlm1)+"\n")
    f.write( "merge time used: " +str(tmerge2-tmerge1)+"\n")
    f.write( "total time used: : "+str(t2-t1))