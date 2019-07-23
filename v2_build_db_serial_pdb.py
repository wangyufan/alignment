import os, sys, time
import pdb2zernike
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx.option_parser import option_parser
from libtbx import easy_pickle
import argparse

banner = "--------------Build Zernike Moment DataBase----------------"


def run(args):
  parser = argparse.ArgumentParser()
  parser.add_argument("--prefix", help="result path", default="myTestDB",type=str)
  parser.add_argument("-path", help="db path",  type=str)
  parser.add_argument("--np", help="number of point covering [0,1]", default=50, type=int)
  parser.add_argument("--fix_dx", help="Whether keeping default dx=0.7A or not",default=True, type=bool)
  parser.add_argument("--nmax", help="nmax", default=20, type=int)
  parser.add_argument("--qmax", help="rmax", default=0.3, type=float)
  args = parser.parse_args()

  path = args.path
  print("filepath:",path)
  nmax = args.nmax
  np = args.np
  fix_dx = args.fix_dx
  prefix = args.prefix
  print("prefix:",prefix)

  pdb_dir = os.listdir(path)
  files = [f for f in pdb_dir if f.endswith('pdb')]
  sorted_files = sorted(files, key=lambda oneFileName:int(oneFileName.split(".")[0]))

  nlm_coefs = []
  codes = []
  for file in sorted_files:
    code = file.split('\n')[0].split('.')[0]
    file = path+file

    mom_obj, vox_obj, pdb = pdb2zernike.zernike_moments(file, nmax=nmax, np=np, fix_dx=fix_dx, coef_out=False, calc_intensity=False )
    if(mom_obj is None):
      print(code, "NOT processed, please check the file")
      continue
    codes.append( code )
    nlm_coefs.append( mom_obj.moments().coefs().deep_copy() )
    print(code, "processed.")

  easy_pickle.dump(prefix+".nlm", nlm_coefs)
  easy_pickle.dump(prefix+".codes", codes)


if __name__ == "__main__":
  t1 = time.time()
  args = sys.argv[1:]
  run(args)
  t2 = time.time()
  print("time used: ", t2-t1)
