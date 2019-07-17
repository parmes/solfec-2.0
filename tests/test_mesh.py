# Solfec-2.0 unit test: ARGV input command
import unittest
import os

output0 = \
'''MESH_0_nodes = [0,0,0,
1,0,0,
1,1,0,
0,1,0,
0,0,1,
1,0,1,
1,1,1,
0,1,1,
0.5,0.5,2,
1,1,3,
0,1,3,
0.5,0.5,3,
0.5,0.5,4]
MESH_0_elements = [8,0,1,2,3,4,5,6,7,0,
5,4,5,6,7,8,0,
6,6,7,8,9,10,11,0,
4,9,10,11,12,0]
MESH_0_colors = [1,4,0,1,2,3,2,
4,4,5,6,7,3,
3,10,11,12,4]
'''

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 mesh.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
