import os, sys
import mol2
from scitbx.math import zernike_align_fft as fft_align
from scitbx.math import *
import time
import model_interface
import multiprocessing
from multiprocessing import Pool
import argparse
import math as xmath
import numpy as np
from libtbx import easy_pickle

def computeCc(process_total, process_n, perSizeLeft, perSizeRight, perSize, listlength, filelist, nmax, rmax):
#compute nlm_cc
  ccArr = []
  start_line_left = process_n*perSizeLeft
  end_line_left = (process_n + 1)*perSizeLeft
  start_line_right = listlength - (process_n+1)*perSizeRight
  end_line_right = listlength - process_n*perSizeRight

  if start_line_right - end_line_left < 0:
    #arrive the stop point, the interval is the last difference of start_line_right and end_line_left
    middle = start_line_right - end_line_left + perSize
    end_line_left = end_line_left - perSizeLeft + int(xmath.ceil(middle / 2))
    start_line_right = end_line_left 

  #if stop before left meet right  
  if((start_line_right - end_line_left > 0) and (process_n == process_total)):
  	end_line_left = start_line_right

  print("process_n, start_left, end_left, start_right, end_right:",process_n, start_line_left, end_line_left, start_line_right, end_line_right)
  for i in (range(start_line_left,end_line_left) + range(start_line_right,end_line_right)):
    for j in range(i+1,listlength): 
      fix_nlm_array = nlm_array(nmax)
      mov_nlm_array = nlm_array(nmax)
      fix_nlm = fix_nlm_array.nlm()
      mov_nlm = mov_nlm_array.nlm()
      fix = filelist[i]
      mov = filelist[j]
      fix_nlm_array.load_coefs(fix_nlm, fix)    
      mov_nlm_array.load_coefs(mov_nlm, mov) 
      # t0 = time.time()
      align_obj = fft_align.align( fix_nlm_array, mov_nlm_array, nmax=nmax,refine=True)
      # t1 = time.time()
      # print("align time used:" , t1-t0)
      cc = align_obj.get_cc()
      # t2 = time.time()
      # print("get cc time used:" , t2-t1)
      ccArr.append((i, j, cc))
      print(i, j, cc)
  return ccArr

def run(args):
  t1 = time.time()
  parser = argparse.ArgumentParser()
  parser.add_argument("-output", help="result path", type=str)
  parser.add_argument("-input", help="input nlm path", type=str)
  parser.add_argument("-processnum", help="processnum number", type=int)
  parser.add_argument("-node", help="node number", type=int)
  parser.add_argument("-node_n", help="node number", type=int)
  parser.add_argument("--nmax", help="nmax", default=20, type=int)
  parser.add_argument("--rmax", help="rmax", type=float)
  args = parser.parse_args()
  nodenum = args.node
  node_n = args.node_n
  nmax = args.nmax
  rmax = args.rmax
  prefix = args.input
  nlm_coefs=easy_pickle.load(prefix+".nlm")
  listlength = len(nlm_coefs)
  targetfile = args.output
  processnum = args.processnum
  perSize = int(0.5+ ( listlength / (nodenum*processnum)))
  perSizeLeft = int(xmath.floor(perSize / 2))
  perSizeRight = perSize - perSizeLeft
  print("perSizeRight, perSizeLeft, perSize, listlength:",perSizeRight, perSizeLeft, perSize, listlength)

  #cc process pool
  tnlm1 = time.time()
  res = []
  result = []
  pool = multiprocessing.Pool(processes = processnum)
  process_left = (node_n-1)*processnum
  process_right = node_n*processnum
  process_total = nodenum*processnum - 1
  print("node_n", node_n, process_left, process_right)
  for process_n in range(process_left, process_right):
    result.append(pool.apply_async(computeCc, (process_total, process_n, perSizeLeft, perSizeRight, perSize, listlength, nlm_coefs, nmax, rmax)))
  pool.close()
  pool.join()
  for item in result:
  	res.extend(item.get())
  np.save(targetfile, res)
  t2 = time.time()
  print("time used", t2 - t1)

if __name__=="__main__":
  args = sys.argv[1:]
  print(args)
  run(args)