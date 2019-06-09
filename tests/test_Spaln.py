#!/usr/bin/env python
# Author: Tomas Bruna

import unittest
import sys
import os
import multiprocessing
import subprocess
import shutil
import filecmp


class TestSpaln(unittest.TestCase):

    def testSpalnToGff(self):
        os.chdir(testDir + "/test_Spaln")
        command = prothint.binDir + "/spaln_to_gff.py < test.spaln --intronScore 10 \
                  --startScore 10 --stopScore 10 --gene 289_g.fasta --prot \
                  156304_0:00222f.fasta > single.gff"
        subprocess.call(command, shell=True)

        self.assertEqual(filecmp.cmp("result_single.gff",
                                     "single.gff"), True)
        os.remove("single.gff")

    def testSpaln(self, pbs=True):
        prothint.workDir = testDir + "/test_Spaln"
        prothint.proteins = prothint.workDir + "/proteins.fasta"
        prothint.runSpaln(prothint.workDir + "/pairs.out", pbs=pbs)

        command = "diff <(sort " + prothint.workDir + "/Spaln/spaln.gff) \
            <(sort " + prothint.workDir + "/result_spaln.gff)"
        diffResult = subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)
        if (diffResult == 0):
            shutil.rmtree(prothint.workDir + "/Spaln")

    def testSpalnNoPbs(self):
        self.testSpaln(pbs=False)


if __name__ == '__main__':
    testDir = os.path.abspath(os.path.dirname(__file__))
    sys.path.append(testDir + "/../bin")

    import prothint
    prothint.threads = str(multiprocessing.cpu_count())
    prothint.binDir = testDir + "/../bin"

    unittest.main()
