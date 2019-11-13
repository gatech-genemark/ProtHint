#!/usr/bin/env python3
# Author: Tomas Bruna

import unittest
import sys
import os
import subprocess


class TestProcessProSplignOutput(unittest.TestCase):

    def compareFiles(self, file1, file2):
        command = "diff <(sort " + file1 + " ) \
                <(sort " + file2 + " )"
        diffResult = subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)

    def testProcessOutput(self):
        prothint.workDir = testDir + "/test_processProSplignOutput"
        prothint.processProSplignOutput()
        os.chdir(prothint.workDir)

        self.compareFiles("ProSplign/prosplign_combined.gff",
                          "test_prosplign_combined.gff")
        self.compareFiles("prothint.gff",
                          "test_prothint.gff")
        self.compareFiles("evidence.gff",
                          "test_evidence.gff")
        self.compareFiles("prothint_augustus.gff",
                          "test_prothint_augustus.gff")
        self.compareFiles("evidence_augustus.gff",
                          "test_evidence_augustus.gff")

        os.remove("ProSplign/prosplign_combined.gff")
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
