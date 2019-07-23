# Solfec-2.0 unit test: FRICTION input command
import unittest, os

output0 = \
'''FRICTION_0_color1 = 1
FRICTION_0_color2 = 2
FRICTION_0_static = 0.2
FRICTION_0_dynamic = 0.3
FRICTION_1_color1 = 0
FRICTION_1_color2 = 0
FRICTION_1_static = 0
FRICTION_1_dynamic = 0
FRICTION_2_color1 = 0
FRICTION_2_color2 = 1
FRICTION_2_static = 0.1
FRICTION_2_dynamic = 0.1
'''

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 friction.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
