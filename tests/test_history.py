# Solfec-2.0 unit test: HISTORY input command
import unittest, os

output0 = \
'''HISTORY_0_entity = 'TIME'
HISTORY_1_entity = '|P|'
HISTORY_1_bodnum = 0
HISTORY_1_points = [(0,0,0)]
HISTORY_2_entity = '|S|'
HISTORY_2_bodnum = 0
HISTORY_2_points = [(1,1,1)]
HISTORY_2_filepath = history_mises.txt
'''

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 history.py')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
