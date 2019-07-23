import sys,os
import multiprocessing
from multiprocessing import Pool
import math

# def read_xyz(xyz_file):
#   xyz_file = '/Users/wyf/Downloads/protein_shape_retrieval_contest/'+ xyz_file + '.off'
#   # xyz_file = '/home/dongxq/Downloads/protein_shape_retrieval_contest/'+ xyz_file + '.off'
#   f = open( xyz_file, 'r')
#   xyz = []
#   for line in f.readlines()[2:]:
#     keys = line.split()
#     if (len(keys) == 3):
#       this_xyz = [ float( keys[0] ), float(keys[1]), float(keys[2]) ]
#       xyz.append( this_xyz )
#   f.close()
#   return xyz
  

# def write_pdb(xyz, filename='converted.pdb'):
#   res_id = 1
#   # filename = '/home/dongxq/Downloads/off2pdb/' + filename
#   filename = '/Users/wyf/Downloads/off2pdb/' + filename
#   out = open( filename, 'w')
#   for x,y,z in xyz:
#     out.write("ATOM  %5d  CA  ASP  %4d    %8.3f%8.3f%8.3f  1.00  1.00\n"%(res_id, res_id, x, y, z))
#     res_id = res_id + 1
#   out.close()

# def get_pdb(process_n, start, end):
#   for i in range(start, end):
#     print(start, end, i)
#     xyz_file = i
#     pdb_file = str(xyz_file) + '.pdb'
#     xyz = read_xyz( xyz_file )
#     write_pdb( xyz, filename=pdb_file )

# if __name__ == "__main__":
#   args = sys.argv[1:]
#   processnum = int(args[0])
#   total = int(args[1])
#   per = int(total / processnum)
#   pool = multiprocessing.Pool(processes = processnum)
#   res = []
#   for process_n in range(processnum):
#     s = per*process_n
#     if(process_n == processnum-1):
#       e = total
#     else:
#       e = per*(process_n+1)
#     # print(process_n, s, e)
#     res.append(pool.apply_async(get_pdb, (process_n, s, e)))
#   pool.close()
#   pool.join()
#   for item in res:
#     print(item)

def read_xyz(xyz_file):
  # xyz_file = '/Users/wyf/Downloads/protein_shape_retrieval_contest/'+ xyz_file + '.off'
  xyz_file = '/home/dongxq/Downloads/protein_shape_retrieval_contest/'+ xyz_file + '.off'
  f = open( xyz_file, 'r')
  xyz = []
  for line in f.readlines()[2:]:
    keys = line.split()
    if (len(keys) == 3):
      this_xyz = [ float( keys[0] ), float(keys[1]), float(keys[2]) ]
      xyz.append( this_xyz )
  f.close()
  return xyz
  

def write_pdb(xyz, filename='converted.pdb'):
  res_id = 1
  filename = '/home/dongxq/Downloads/off2pdb/' + filename
  # filename = '/Users/wyf/Downloads/off2pdb/' + filename
  out = open( filename, 'w')
  for x,y,z in xyz:
    out.write("ATOM  %5d  CA  ASP  %4d    %8.3f%8.3f%8.3f  1.00  1.00\n"%(res_id, res_id, x, y, z))
    res_id = res_id + 1
  out.close()


if __name__ == "__main__":
  args = sys.argv[1:]
  xyz_file = args[0]
  pdb_file = xyz_file +'.pdb'
  xyz = read_xyz( xyz_file )
  write_pdb( xyz, filename=pdb_file )
  # if( len(args) == 1 ):
  #   xyz_file = args[0]
  #   pdb_file = xyz_file.split('.')[0]+'.pdb'
  # else:
  #   xyz_file = args[0]
  #   pdb_file = args[1]


