# Solfec-2.0 unit test: RESTAIN input command
import unittest, os

output0 = \
'''RESTRAIN_0_bodnum = 0
RESTRAIN_0_points = [(1,0,0)]
RESTRAIN_0_color = 0
RESTRAIN_1_bodnum = 0
RESTRAIN_1_points = [(0,0,0),(1,0,0)]
RESTRAIN_1_color = 0
RESTRAIN_2_bodnum = 0
RESTRAIN_2_points = 
RESTRAIN_2_color = 4
'''

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 restrain.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
