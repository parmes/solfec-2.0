# Solfec-2.0 input parsing test: insertions, runs, and deletions
import sys
sys.path.append('../python')
from mesh_hex import MESH_HEX

n = 4
m = 4
matnum = MATERIAL (1E3, 1E9, 0.25, 0.0)
bodnum_list = []

def insert_meshes(n, m):
  for i in range(0,n):
    for j in range(0,n):
      for k in range(0,n):
        nodes = [(0+i, 0+j, 0+k),
                 (1+i, 0+j, 0+j),
                 (1+i, 1+j, 0+k),
                 (0+i, 1+j, 0+k),
                 (0+i, 0+j, 1+k),
                 (1+i, 0+j, 1+k),
                 (1+i, 1+j, 1+k),
                 (0+i, 1+j, 1+k)]
        colors = [1, 2, 3, 4, 5, 6, 7, 8]
        bodnum = MESH_HEX (nodes, m, m, m, matnum, colors)
        bodnum_list.append(bodnum)

if len(ARGV()) == 0:
  print('INFO: neither INS nor DEL arguments passed')
else:
  for arg in ARGV():
    if arg == 'INS':
      insert_meshes(n, m)
      print('INSERTION')
      for num in bodnum_list:
        print_MESH(num)
      RUN (1.0, 1.0)
