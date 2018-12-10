import os, sys, time
import pdb2zernike
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx.option_parser import option_parser
from libtbx import easy_pickle
import argparse

banner = "--------------Build Zernike Moment DataBase----------------"

def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "\n  Usage: libtbx.python build_db.py path=path nmax=nmax np=np_of_point fix_dx=True/False"
  print >> out, "\n  Attention: \n  The path should be a directory that contains ONLY PDB files\n\n"


def read(path='./'):
  files = os.listdir(path)
  return files

def run(args):

  # parser = argparse.ArgumentParser()
  # parser.add_argument("--prefix", help="result path", default="myTestDB",type=str)
  # parser.add_argument("-path", help="db path",  type=str)
  # parser.add_argument("--np", help="number of point covering [0,1]", default=50, type=int)
  # parser.add_argument("--fix_dx", help="Whether keeping default dx=0.7A or not",default=True, type=bool)
  # parser.add_argument("--nmax", help="nmax", default=20, type=int)
  # parser.add_argument("--qmax", help="rmax", default=0.3, type=float)
  # args = parser.parse_args()


  path = args.path
  print "filepath:",path
  nmax = args.nmax
  np = args.np
  fix_dx = args.fix_dx
  prefix = args.prefix
  print "prefix:",prefix


  files = read(path)
  nlm_coefs = []
  nl_coefs = []
  nn_coefs = []
  codes = []
  rmax = []
  for file in files[1:]:
    code = file.split('\n')[0].split('.')[0]
    file = path+file

    mom_obj, vox_obj, pdb = pdb2zernike.zernike_moments(file, nmax=nmax, np=np, fix_dx=fix_dx, coef_out=False, calc_intensity=False )
    if(mom_obj is None):
      print code, "NOT processed, please check the file"
      continue
    codes.append( code )
    rmax.append( vox_obj.rmax() )
    nlm_coefs.append( mom_obj.moments().coefs().deep_copy() )
    nn_coefs.append( mom_obj.fnn().coefs().deep_copy() )
    nl_coefs.append( mom_obj.fnl().coefs() )
    print code, "processed."

  easy_pickle.dump(prefix+".nlm", nlm_coefs)
  easy_pickle.dump(prefix+".nn", nn_coefs)
  easy_pickle.dump(prefix+".rmax", rmax)
  easy_pickle.dump(prefix+".codes", codes)
  easy_pickle.dump(prefix+".nl", nl_coefs)


if __name__ == "__main__":
  t1 = time.time()
  args = sys.argv[1:]
  run(args)
  t2 = time.time()
  print "time used: ", t2-t1
