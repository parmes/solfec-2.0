# Solfec-2.0 unit test: MESH input command
import unittest, os, random, glob, ast, itertools, bisect

# https://stackoverflow.com/questions/38346013/binary-search-in-a-python-list
def search(alist, item):
    'Locate the leftmost value exactly equal to item'
    i = bisect.bisect_left(alist, item)
    if i != len(alist) and alist[i] == item:
        return i
    raise ValueError

def read_stdout_meshes(output):
  stdout_meshes = dict()
  nmsh = output.count('_nodes')
  jn = 0
  je = 0
  for i in range(0,nmsh):
    jn0 = jn
    jn = output.find('_nodes', jn+1)
    j0 = output.rfind('_',jn0, jn-1)
    bodnum = int(output[j0+1:jn])
    k = output.find('[', jn)
    l = output.find(']', jn)
    nlist = ast.literal_eval(output[k:l+1])
    nlist0 = [tuple(nlist[i:i+3]) for i in range(0, len(nlist), 3)]
    nlist = sorted(nlist0)
    je = output.find('_elements', je+1)
    k = output.find('[', je)
    l = output.find(']', je)
    elist = ast.literal_eval(output[k:l+1])
    el = len(elist)
    ed = [0]
    k = 0
    while k < el:
      l = elist[k]+1
      elist[k+1:k+l] = [search(nlist,nlist0[i]) for i in elist[k+1:k+l]]
      k += l+1
      ed.append(k)
    elist = [tuple(elist[ed[i]:ed[i+1]]) for i in range(0, len(ed)-1)]
    elist = sorted(elist)
    stdout_meshes[bodnum] = (nlist, elist) # output sorted nodes and elements
  return stdout_meshes
  
def read_debug_meshes(np):
  debug_meshes = dict()
  rank_nodes = []
  body_elems = dict()
  for r in range(0,np): # read nodes
    with open('debug_nodes%d.txt'%r) as f:
      lines = [line.rstrip('\n') for line in f]
    assert 'COUNT' in lines[1]
    nnod = int(lines[1].split()[1])
    nlist = [tuple(float(f) for f in l.split()) for l in lines[2:2+nnod]]
    rank_nodes.append(nlist)
  for r in range(0,np): # read elements
    with open('debug_elements%d.txt'%r) as f:
      lines = [line.rstrip('\n') for line in f]
    assert 'COUNT' in lines[1]
    nele = int(lines[1].split()[1])
    elist = []
    i = 2
    j = 0
    while j < nele:
      assert 'BODNUM' in lines[i]
      bodnum = int(lines[i].split()[1])
      assert 'MATNUM' in lines[i+1]
      matnum = int(lines[i+1].split()[1])
      assert 'TYPE' in lines[i+2]
      etype = int(lines[i+2].split()[1])
      eledef = [etype]
      i += 3
      for k in range(0,etype):
        tmp = lines[i].split()
        rng = int(tmp[0])
        idx = int(tmp[1])
        i += 1
        eledef.append((rng, idx))
      eledef.append(matnum) # (type, (r0,i0), ... (rtype-1,itype-1), matnum)
      body_elems.setdefault(bodnum, []).append(eledef)
      j += 1
  for bodnum in body_elems: # process elements
    i = 0
    nmap = dict() # node (r,i) to linear index
    for ele in body_elems[bodnum]:
      for nd in ele[1:1+ele[0]]:
        if not nd in nmap:
          nmap[nd] = i
          i += 1
    nlist = [0]*i # final node list
    for nd in nmap:
      nlist[nmap[nd]] = rank_nodes[nd[0]][nd[1]]
    nlist = sorted(nlist)
    elist = [] # final element list
    for ele in body_elems[bodnum]:
      eledef = [ele[0]]
      for nd in ele[1:1+ele[0]]:
        pt = rank_nodes[nd[0]][nd[1]]
        eledef.append(search(nlist, pt))
      eledef.append(ele[-1])
      elist.append(tuple(eledef))
    elist = sorted(elist)
    debug_meshes[bodnum] = (nlist, elist) # output sorted nodes and elements
  return debug_meshes

class CompareMeshesAssertions:
  def assertSameMeshes(self, stdout, debug):
    if len(stdout) != len(debug):
      raise AssertionError('Mesh counts differ: %d != %d' % (len(stdout), len(debug)))

    for bodnum in stdout:
      if not bodnum in debug:
        raise AssertionError('Mesh bodnum %d missing in debug files' % bodnum)

    for bodnum in debug:
      if not bodnum in stdout:
        raise AssertionError('Mesh bodnum %d missing in stdout' % bodnum)

    for bodnum in stdout:
      nodes0 = stdout[bodnum][0]
      nodes1 = debug[bodnum][0]
      if len(nodes0) != len(nodes1):
        raise AssertionError('For bodnum %d node counts differ: %d (stdout) != %d (debug)' % (bodnum, len(nodes0), len(nodes1)))

      for i in range(0, len(nodes0)):
        n0 = nodes0[i]
        n1 = nodes1[i]
        if n0 != n1:
          raise AssertionError('For bodnum %d node %d differs: %s (stdout) != %s (debug)' % (bodnum, i, n0, n1))

      eles0 = stdout[bodnum][1]
      eles1 = debug[bodnum][1]
      if len(eles0) != len(eles1):
        raise AssertionError('For bodnum %d element counts differ: %d (stdout) != %d (debug)' % (bodnum, len(eles0), len(eles1)))

      for i in range(0, len(eles0)):
        e0 = eles0[i]
        e1 = eles1[i]
        if e0 != e1:
          raise AssertionError('For bodnum %d node %d differs: %s (stdout) != %s (debug)' % (bodnum, i, e0, e1))

def run_test(arg_string):
  np = random.randint(2,12)
  cmd = 'mpirun -np %d --oversubscribe ../solfec4 --debug_print insertions_and_deletions.py ' % np
  solfec = os.popen(cmd + arg_string)
  output = solfec.read()
  solfec.close()

  stdout_meshes = read_stdout_meshes(output)
  debug_meshes = read_debug_meshes(np)

  filelist = glob.glob('debug_*.txt')
  for filepath in filelist:
    try:
      os.remove(filepath)
    except:
      print("Error while deleting file : ", filepath)

  return (stdout_meshes, debug_meshes)

class test_MESH_i(unittest.TestCase, CompareMeshesAssertions):
  def test_MESH_i(self):
    print('\ntesting MESH insertions')
    stdout_meshes, debug_meshes = run_test('INSMESH')
    self.assertSameMeshes(stdout_meshes, debug_meshes)

class test_MESH_id(unittest.TestCase, CompareMeshesAssertions):
  def test_MESH_id(self):
    print('\ntesting MESH insertions-deletions')
    stdout_meshes, debug_meshes = run_test('INSMESH DEL')
    self.assertSameMeshes(stdout_meshes, debug_meshes)

class test_MESH_idi(unittest.TestCase, CompareMeshesAssertions):
  def test_MESH_idi(self):
    print('\ntesting MESH insertions-deletions-insertions')
    stdout_meshes, debug_meshes = run_test('INSMESH DEL INSMESH')
    self.assertSameMeshes(stdout_meshes, debug_meshes)

class test_MESH_idid(unittest.TestCase, CompareMeshesAssertions):
  def test_MESH_idid(self):
    print('\ntesting MESH insertions-deletions-insertions-deletions')
    stdout_meshes, debug_meshes = run_test('INSMESH DEL INSMESH DEL')
    self.assertSameMeshes(stdout_meshes, debug_meshes)

class test_MESH_ididi(unittest.TestCase, CompareMeshesAssertions):
  def test_MESH_ididi(self):
    print('\ntesting MESH insertions-deletions-insertions-deletions-insertions')
    stdout_meshes, debug_meshes = run_test('INSMESH DEL INSMESH DEL INSMESH')
    self.assertSameMeshes(stdout_meshes, debug_meshes)
