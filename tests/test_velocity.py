# Solfec-2.0 unit test: VELOCITY input command
import unittest, os

output0 = \
'''VELOCITY_0_bodnum = 0
VELOCITY_0_linear = (1,0,0)
VELOCITY_0_angular = (0,0,0)
VELOCITY_1_bodnum = 0
VELOCITY_1_linear = (0,0,0)
VELOCITY_1_angular = (1,0,0)
'''

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 velocity.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
