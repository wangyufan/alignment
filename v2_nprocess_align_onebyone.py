import os, sys
import mol2
from scitbx.math import zernike_align_fft as fft_align
import time
import model_interface
import multiprocessing
from multiprocessing import Pool
import argparse
import math as xmath
import numpy as np

def computeCc(process_n, perSizeLeft, perSizeRight, perSize, listlength, filelist, nmax, rmax):
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

  print "process_n, start_left, end_left, start_right, end_right:",process_n, start_line_left, end_line_left, start_line_right, end_line_right
  for i in (range(start_line_left,end_line_left) + range(start_line_right,end_line_right)):
    for j in range(i+1,listlength):      
      fix_model = model_interface.build_model( filelist[i], "pdb", nmax, rmax )
      mov_model = model_interface.build_model( filelist[j], "pdb", nmax, rmax )
      fix_nlm_array = fix_model.nlm_array
      mov_nlm_array = mov_model.nlm_array
      align_obj = fft_align.align( fix_nlm_array, mov_nlm_array, nmax=nmax,refine=True)
      cc = align_obj.get_cc()
      ccArr.append((i, j, cc))
  return ccArr

def run(args):
  t1 = time.time()
  parser = argparse.ArgumentParser()
  parser.add_argument("-output", help="result path", type=str)
  parser.add_argument("-input", help="input mol2/pdb path", type=str)
  parser.add_argument("-processnum", help="processnum number", type=int)
  parser.add_argument("--nmax", help="nmax", default=20, type=int)
  parser.add_argument("--rmax", help="rmax", type=float)
  args = parser.parse_args()
  nmax = args.nmax
  rmax = args.rmax
  filepath = args.input
  filelist = []
  # sorted filedir
  files = os.listdir(filepath)
  pdbs = [f for f in files if f.endswith('pdb')]
  filedir = sorted(pdbs, key=lambda oneFileName:int(oneFileName.split(".")[0]))
  print filedir
  for filename in filedir:
    filelist.append(os.path.join('%s%s' % (filepath, filename)))
  listlength = len(filelist)
  targetfile = args.output
  processnum = args.processnum
  perSize = int(xmath.ceil( listlength / processnum))
  perSizeLeft = int(xmath.ceil(perSize / 2))
  perSizeRight = perSize - perSizeLeft
  print "perSizeRight:",perSizeRight, perSizeLeft, perSize, listlength

  #cc process pool
  tnlm1 = time.time()
  res = []
  result = []
  pool = multiprocessing.Pool(processes = processnum)
  for process_n in range(processnum):
    result.append(pool.apply_async(computeCc, (process_n, perSizeLeft, perSizeRight, perSize, listlength, filelist, nmax, rmax)))
  pool.close()
  pool.join()
  for item in result:
    res.extend(item.get())
  np.save(targetfile, res)
  # with open(targetfile,"w") as f:
  #   for line in res:
  #     f.write(str(line)+'\n')
  # f.close()
  t2 = time.time()
  print "time used", t2 - t1

if __name__=="__main__":
  args = sys.argv[1:]
  print args
  run(args)