# Solfec-2.0 unit test: ARGV input command
import unittest
import os

class test(unittest.TestCase):
  def test(self):
    solfec = os.popen('../solfec4 argv.py Non Solfec Args')
    output = solfec.read()
    solfec.close()
    self.assertEqual(output, "['Non', 'Solfec', 'Args']\n")
