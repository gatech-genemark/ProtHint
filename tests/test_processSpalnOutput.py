#!/usr/bin/env python
# Author: Tomas Bruna

import unittest
import sys
import os
import subprocess
import shutil

class TestProcessSpalnOutput(unittest.TestCase):

    def compareFiles(self, file1, file2):
        command = "diff <(sort " + file1 + " ) \
                <(sort " + file2 + " )"
        diffResult = subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)

    def testProcessOutput(self):
        prothint.workDir = testDir + "/test_processSpalnOutput"
        os.chdir(prothint.workDir)

        shutil.copyfile("Spaln/spaln.gff", "Spaln/spaln_orig.gff")
        prothint.processSpalnOutput("diamond/diamond.out")
        shutil.move("Spaln/spaln_orig.gff", "Spaln/spaln.gff")

        self.compareFiles("prothint.gff",
                          "test_prothint.gff")
        self.compareFiles("evidence.gff",
                          "test_evidence.gff")
        self.compareFiles("prothint_augustus.gff",
                          "test_prothint_augustus.gff")
        self.compareFiles("evidence_augustus.gff",
                          "test_evidence_augustus.gff")

        os.remove("prothint.gff")
        os.remove("evidence.gff")
        os.remove("prothint_augustus.gff")
        os.remove("evidence_augustus.gff")

if __name__ == '__main__':
    testDir = os.path.abspath(os.path.dirname(__file__))
    sys.path.append(testDir + "/../bin")

    import prothint
    prothint.binDir = testDir + "/../bin"

    unittest.main()
