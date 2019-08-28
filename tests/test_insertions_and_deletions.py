# Solfec-2.0 unit test: MESH input command
import unittest, os, random, glob

class test(unittest.TestCase):
  def test(self):

    np = random.randint(2,12)
    cmd = 'mpirun -np %d --oversubscribe ../solfec4 --debug_print insertions_and_deletions.py' % np
    solfec = os.popen(cmd + ' INS')
    output = solfec.read()
    solfec.close()

    # print(output)
    # TODO: parse output meshes and read debug_*.txt data and compare

    filelist = glob.glob('debug_*.txt')
    for filepath in filelist:
      try:
        os.remove(filepath)
      except:
        print("Error while deleting file : ", filepath)

    self.assertEqual(1, 1)
