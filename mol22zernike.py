import sys, os
import math as smath
from scitbx.array_family import flex
import iotbx.phil
from iotbx import pdb, ccp4_map
import time
from mmtbx.monomer_library import server, pdb_interpretation

from scitbx import math
from sastbx import zernike_model as zm
from iotbx.xplor import map as xplor_map
from cctbx import uctbx, sgtbx

from libtbx import easy_pickle
from sastbx.interface import get_input
import mol2


master_params = iotbx.phil.parse("""\
zernike{
  mol2_file = None
  .type=path
  .help="the mol2 file"

  qmax = 0.3
  .type=float
  .help="maximum q value, for which the intensity to be evaluated"

  nmax=20
  .type=int
  .help="maximum order of zernike expansion"

  np=50
  .type=int
  .help="number of point covering [0,1]"

  uniform = True
  .type=bool
  .help="Whether the uniform density is used"

  fix_dx = True
  .type=bool
  .help="Whether the dx will be adjusted to match the size of molecule"

  buildmap=True
  .type=bool
  .help="Whether xplor map will be constructed or not"

  shift = False
  .type=bool
  .help="Whether the pdb coordates will be shifted or not"

  coef_out=True
  .type=bool
  .help="Whether dump zernike moments to picle files"
}
""")


banner = "--------------Zernike Moments Calculation and Map Construction----------------"

def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "Usage: libtbx.python pdb2zernike.py pdbfile=pdbfile nmax=nmax buildmap=True/False shift=True/False np=np_of_point"

def ccp4_map_type(map, N, radius,file_name='map.ccp4'):
  grid = flex.grid(N*2+1, N*2+1,N*2+1)
  map.reshape( grid )
  ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=uctbx.unit_cell(" %s"%(radius*2.0)*3+"90 90 90"),
      space_group=sgtbx.space_group_info("P1").group(),
      gridding_first=(0,0,0),
      gridding_last=(N*2, N*2, N*2),
      map_data=map,
      labels=flex.std_string(["generated from zernike moments"]))



def xplor_map_type(m,N,radius,file_name='map.xplor'):
  gridding = xplor_map.gridding( [N*2+1]*3, [0]*3, [2*N]*3)
  grid = flex.grid(N*2+1, N*2+1,N*2+1)
  m.reshape( grid )
  uc = uctbx.unit_cell(" %s"%(radius*2.0)*3+"90 90 90")
  xplor_map.writer( file_name, ['no title lines'],uc, gridding,m) # is_p1_cell=True)  # True)


def zernike_moments(mol2_file, nmax=20, np=50, uniform=True, fix_dx=False, np_on_grid=20, shift=False, buildmap=False, coef_out=True, calc_intensity=True, external_rmax=-1):
  base = mol2_file.split('.')[0]
  splat_range = 0
  fraction = 0.9
  default_dx = 0.5

  models = mol2.read_Mol2_file( mol2_file )
  if(len( models )== 0):
    return None,None,None
  ## for testing purpose, here only the first model in the mol2 file is converted ##
  this_molecule = models[0]
  ## models should be a list of molecules, same ligand at different conformations ##

  all_atoms = this_molecule.atom_list
  if(len( all_atoms )== 0):
    return None,None,None

  # predefine some arrays we will need
  xyz = flex.vec3_double()
  charges = flex.double()
  # keep track of the atom types we have encountered
  for atom in all_atoms:
    xyz.append( ( atom.X, atom.Y, atom.Z ) )
    charges.append( atom.Q )

  if(xyz.size() == 0):
    return None,None,None
  if (uniform):
    density=flex.double(xyz.size(),1.0)
  else: # use charges
    density=charges

  voxel_obj = math.sphere_voxel(np,splat_range,uniform,fix_dx,external_rmax, default_dx, fraction,xyz,density)
  np = voxel_obj.np()
  rmax=voxel_obj.rmax()/fraction

  #print base, "RMAX: ", voxel_obj.rmax()

  #print time.time()
  grid_obj = math.sphere_grid(np, nmax)
  pdb_out = False
  grid_obj.clean_space(voxel_obj, pdb_out)
  #print time.time(), "END GRID CONST"
  grid_obj.construct_space_sum()

  #print time.time(), "SS CALC"
  mom_obj = math.zernike_moments(grid_obj, nmax)
  #print time.time(), "END MOM CALC"

  moments = mom_obj.moments()
  if(coef_out):
    nn = mom_obj.fnn()
    easy_pickle.dump(base+'.nlm.pickle', moments.coefs() )
    easy_pickle.dump(base+'.nn.pickle', nn.coefs() )

  #
####### The following will be optional ###
  if(shift):
    shift = [rmax, rmax, rmax]
    centered_xyz = voxel_obj.xyz() + shift
    out_mol2_filename =base+'_centered.mol2'
    for a,xyz in zip( all_atoms, centered_xyz):
      a.X = xyz[0]
      a.Y = xyz[1]
      a.Z = xyz[2]

    mol2.write_mol2( this_molecule, out_mol2_filename )

    original_map = voxel_obj.map()
    xplor_map_type( original_map, np, rmax, file_name=base+'_pdb.xplor')


######## Calculate Intnesity Profile ############
  if(calc_intensity):
    q_array = flex.double( range(51) )/100.0
    z_model = zm.zernike_model( moments, q_array, rmax, nmax)
    nn = mom_obj.fnn()
    intensity = z_model.calc_intensity(nn)
    #iq_file = open(base+"_"+str(nmax)+"_"+str(int(rmax*fraction))+".zi", 'w')
    iq_file = open(base+"_"+str(nmax)+".zi", 'w')
    for qq, ii in zip(q_array, intensity):
      print>>iq_file, qq,ii

    iq_file.close()
####### END of Intensity Calculation ###########

  if(buildmap):
    print "grid setting up..."
    zga = math.zernike_grid(np_on_grid, nmax, False)
    print "grid setting up...Done"
    zga.load_coefs(moments.nlm(), moments.coefs() )
    print "start reconstruction"
    map=flex.abs(zga.f() )
    print "finished reconstruction"
    xplor_map_type( map, np_on_grid, rmax, file_name=base+'.xplor')
    ccp4_map_type( map, np_on_grid, rmax, file_name=base+'.ccp4' )

  return mom_obj, voxel_obj, models

def disp_inv(inva):
  for indx, inv in zip(inva.nl(), inva.coefs() ):
    print indx, inv

def display(moments):
  for indx, mom in zip(moments.nlm(), moments.coefs()):
    print indx, mom

def run(args):
  params = get_input(args, master_params, "zernike", banner, help)
  if params is None:
    return
  mol2_file = params.zernike.mol2_file
  nmax=params.zernike.nmax
  np = params.zernike.np
  uniform = params.zernike.uniform
  fix_dx = params.zernike.fix_dx
  shift = params.zernike.shift
  buildmap = params.zernike.buildmap
  coef_out = params.zernike.coef_out

  zernike_moments(mol2_file, nmax, np, uniform=uniform, fix_dx=fix_dx, shift=shift, buildmap=buildmap, coef_out=coef_out, calc_intensity=True)


if __name__ == "__main__":
   args=sys.argv[1:]
   t0 = time.time()
   run(args)
   t1 = time.time()
   print "time used %e seconds"%(t1-t0)
