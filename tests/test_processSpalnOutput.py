#!/usr/bin/env python
# Author: Tomas Bruna

import unittest
import sys
import os
import subprocess
import shutil
import common

class TestProcessSpalnOutput(unittest.TestCase):

    def testProcessOutput(self):
        prothint.workDir = testDir + "/test_processSpalnOutput"
        os.chdir(prothint.workDir)

        shutil.copyfile("Spaln/spaln.gff", "Spaln/spaln_orig.gff")
        prothint.processSpalnOutput("diamond/diamond.out")
        shutil.move("Spaln/spaln_orig.gff", "Spaln/spaln.gff")

        self.assertEqual(common.compareFiles("prothint.gff", "test_prothint.gff"), 0)
        self.assertEqual(common.compareFiles("evidence.gff", "test_evidence.gff"), 0)
        self.assertEqual(common.compareFiles("prothint_augustus.gff", "test_prothint_augustus.gff"), 0)
        self.assertEqual(common.compareFiles("evidence_augustus.gff", "test_evidence_augustus.gff"), 0)

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
