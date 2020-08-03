#!/usr/bin/env python3
# Author: Tomas Bruna

import unittest
import sys
import os
import multiprocessing
import subprocess
import shutil
import filecmp


class TestSpaln(unittest.TestCase):

    def testSpaln(self, pbs=True):
        prothint.workDir = testDir + "/test_Spaln"
        # Proteins file is deleted at the end of prothint
        shutil.copy(prothint.workDir + "/proteins.fasta", prothint.workDir + "/proteins_copy.fasta")
        prothint.proteins = prothint.workDir + "/proteins_copy.fasta"
        prothint.runSpaln(prothint.workDir + "/pairs.out", pbs, 25)

        command = "diff <(sort " + prothint.workDir + "/Spaln/spaln.gff) \
            <(sort " + prothint.workDir + "/result_spaln.gff)"
        diffResult = subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)
        if (diffResult == 0):
            shutil.rmtree(prothint.workDir + "/Spaln")

    def testSpalnNoPbs(self):
        self.testSpaln(pbs=False)

    def testSpaln1000(self):
        prothint.workDir = testDir + "/test_Spaln"
        # Proteins file is deleted at the end of prothint
        shutil.copy("prothint.workDir + /proteins.fasta", "prothint.workDir + /proteins_copy.fasta")
        prothint.proteins = prothint.workDir + "/proteins_copy.fasta"
        prothint.runSpaln(prothint.workDir + "/pairs.out", True, 1000)

        command = "diff <(sort " + prothint.workDir + "/Spaln/spaln.gff) \
            <(sort " + prothint.workDir + "/result_spaln_1000.gff)"
        diffResult = subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)
        if (diffResult == 0):
            shutil.rmtree(prothint.workDir + "/Spaln")


if __name__ == '__main__':
    testDir = os.path.abspath(os.path.dirname(__file__))
    sys.path.append(testDir + "/../bin")

    import prothint
    prothint.threads = str(multiprocessing.cpu_count())
    prothint.binDir = testDir + "/../bin"

    unittest.main()
