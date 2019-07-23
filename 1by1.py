import os, sys
import mol2
from scitbx.math import zernike_align_fft as fft_align
import time
import model_interface
import multiprocessing
from multiprocessing import Pool
import argparse
import math as xmath

def computeCc(filelist, nmax, rmax):
#compute nlm_cc    
  fix_model = model_interface.build_model( filelist[0], "pdb", nmax, rmax )
  mov_model = model_interface.build_model( filelist[1], "pdb", nmax, rmax )
  fix_nlm_array = fix_model.nlm_array
  mov_nlm_array = mov_model.nlm_array
  align_obj = fft_align.align( fix_nlm_array, mov_nlm_array, nmax=nmax,refine=True)
  cc = align_obj.get_cc()
  print("when nmax:"+str(nmax)+"----------- cc:"+str(cc))
  return cc

def run(args):
  t1 = time.time()
  parser = argparse.ArgumentParser()
  # parser.add_argument("-output", help="result path", type=str)
  parser.add_argument("-input", help="input mol2/pdb path", type=str)
  parser.add_argument("--nmax", help="nmax", default=20, type=int)
  parser.add_argument("--rmax", help="rmax", type=float)
  args = parser.parse_args()
  nmax = args.nmax
  rmax = args.rmax
  filepath = args.input
  filelist = []
  filedir = os.listdir(filepath)
  print filedir
  for filename in filedir:
    filelist.append(os.path.join('%s%s' % (filepath, filename)))
  data = computeCc(filelist, nmax, rmax)
  return data

if __name__=="__main__":
  args = sys.argv[1:]
  run(args)