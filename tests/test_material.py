# Solfec-2.0 unit test: MATERIAL input command
import unittest, os

output0 = \
'''MATERIAL_0_density = 1000
MATERIAL_0_young = 1e+09
MATERIAL_0_poisson = 0.25
MATERIAL_0_viscosity = 0
'''

class test(unittest.TestCase):
  def test(self):
    print('\ntesting MATERIAL command')
    solfec = os.popen('../solfec4 material.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
