# Solfec-2.0 unit test: PRESCRIBE input command
import unittest, os

output0 = \
'''PRESCRIBE_0_bodnum = 0
PRESCRIBE_0_points = [(1,0,0)]
PRESCRIBE_0_color = 0
PRESCRIBE_0_linear = callback
PRESCRIBE_1_bodnum = 0
PRESCRIBE_1_points = [(0,0,0),(1,0,0)]
PRESCRIBE_1_color = 0
PRESCRIBE_1_linear = callback
PRESCRIBE_2_bodnum = 0
PRESCRIBE_2_points = 
PRESCRIBE_2_color = 4
PRESCRIBE_2_angular = callback
PRESCRIBE_3_bodnum = 0
PRESCRIBE_3_points = 
PRESCRIBE_3_color = 4
PRESCRIBE_3_angular = callback
'''

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 prescribe.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
