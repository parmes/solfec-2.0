# Solfec-2.0 unit test: OUTPUT input command
import unittest, os

output0 = \
'''OUTPUT_0_entities = ['AREA','BPAIR','CF','CFN','CFT','COLOR','CPAIR','DISPL','LINVEL','NUMBER','STRESS']
OUTPUT_0_modes = ['CD','EL','MESH','SURF']
OUTPUT_1_entities = ['COLOR','DISPL']
OUTPUT_1_subset = [0]
OUTPUT_1_modes = ['CD','EL','MESH','SURF']
OUTPUT_2_entities = ['LINVEL','STRESS']
OUTPUT_2_subset = [1]
OUTPUT_2_modes = ['MESH']
OUTPUT_3_entities = ['AREA','CF','CPAIR']
OUTPUT_3_subset = [0,1]
OUTPUT_3_modes = ['CD']
'''

class test(unittest.TestCase):
  def test(self):
    print('\ntesting OUTPUT command')
    solfec = os.popen('../solfec4 output.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
