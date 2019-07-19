# Solfec-2.0 unit test: ARGV input command
import unittest, os

output0 = \
'''['Non', 'Solfec', 'Args']
'''

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 argv.py Non Solfec Args')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, output0)
