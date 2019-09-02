# Solfec-2.0 unit test: GRAVITY input command
import unittest, os

output0 = \
'''GRAVITY_gx = 0.000000e+00
GRAVITY_gy = callback
GRAVITY_gz = 0
'''

class test(unittest.TestCase):
  def test(self):
    print('\ntesting GRAVITY command')
    solfec = os.popen('../solfec4 gravity.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
