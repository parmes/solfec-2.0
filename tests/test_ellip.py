# Solfec-2.0 unit test: ELLIP input command
import unittest, os

output0 = \
'''ELLIP_0_material = 0
ELLIP_0_center = (0,0,0)
ELLIP_0_radius = (1,2,3)
ELLIP_0_color = 1
ELLIP_1_material = 0
ELLIP_1_center = (1,1,1)
ELLIP_1_radius = (3,2,1)
ELLIP_1_color = 2
'''

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 ellip.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
