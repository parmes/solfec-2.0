# Solfec-2.0 input parsing test: insertions, runs, and deletions
import sys, random
sys.path.append('../python')
from mesh_hex import MESH_HEX
from static_vars import static_vars

matnum = MATERIAL (1E3, 1E9, 0.25, 0.0)
bodnum_list = []

@static_vars(shift=0)
def insert_meshes():
  n = random.randint(2, 4)
  m = random.randint(2, 4)
  s = insert_meshes.shift # coordinate uniqueness
  insert_meshes.shift += 2*n+1 # separation between mesh blocks
  for i in range(s,s+2*n,2): # individual meshes separated by 1.0
    for j in range(s,s+2*n,2):
      for k in range(s,s+2*n,2):
        nodes = [(0+i, 0+j, 0+k),
                 (1+i, 0+j, 0+k),
                 (1+i, 1+j, 0+k),
                 (0+i, 1+j, 0+k),
                 (0+i, 0+j, 1+k),
                 (1+i, 0+j, 1+k),
                 (1+i, 1+j, 1+k),
                 (0+i, 1+j, 1+k)]
        colors = [1, 2, 3, 4, 5, 6, 7, 8]
        bodnum = MESH_HEX (nodes, m, m, m, matnum, colors)
        bodnum_list.append((bodnum, 'MESH'))

@static_vars(shift=0)
def insert_ellips():
  n = random.randint(4, 8)
  s = insert_ellips.shift # coordinate uniqueness
  insert_ellips.shift += n # shift ellipsoids blocks
  for i in range(s,s+n):
    for j in range(s,s+n):
      for k in range(s,s+n):
        dx = 0.4+random.uniform(0.0, 0.2)
        dy = 0.4+random.uniform(0.0, 0.2)
        dz = 0.4+random.uniform(0.0, 0.2)
        center = (i+dx, j+dy, k+dz)
        a = random.uniform(0.2, 0.4)
        b = random.uniform(0.2, 0.4)
        c = random.uniform(0.2, 0.4)
        radius = (a, b, c)
        color = 1
        bodnum = ELLIP(center, radius, matnum, color)
        bodnum_list.append((bodnum, 'ELLIP'))

def delete_bodies():
  global bodnum_list
  if len(bodnum_list) > 0:
    n = random.randint(1, len(bodnum_list))
    idx = random.sample(range(0, len(bodnum_list)), n)
    todel = [bodnum_list[i] for i in idx]
    for item in todel: DELETE (item[0], item[1])
    bodnum_list = [bodnum_list[i] for i in range(len(bodnum_list)) if i not in idx]

if len(ARGV()) == 0:
  print('INFO: neither of {INSMESH, INSELL, DEL} arguments passed')
else:
  for arg in ARGV():
    if arg == 'INSMESH':
      insert_meshes()
      RUN (0.0, 1.0)
    if arg == 'INSELL':
      insert_ellips()
      RUN (0.0, 1.0)
    if arg == 'DEL':
      delete_bodies()
      RUN (0.0, 1.0)

  # print the current state
  '''
  for item in bodnum_list:
    if item[1] == 'MESH': print_MESH(item[0])
    else: print_ELLIP(item[0])
  '''
