#!/usr/bin/env python
# Author: Tomas Bruna

import unittest
import sys
import os
import multiprocessing
import subprocess
import shutil


class TestProSplign(unittest.TestCase):

    def testProSplign(self):
        prothint.workDir = testDir + "/test_ProSplign/pbs"
        prothint.proteins = prothint.workDir + "/proteins.fasta"
        prothint.runProSplign(pbs=True)

        command = "diff <(sort " + prothint.workDir + "/ProSplign/scored_introns.gff) \
            <(sort " + prothint.workDir + "/result_scored_introns.gff)"
        diffResult = subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)

        command = "diff <(sort " + prothint.workDir + "/ProSplign/prosplign.gff) \
            <(sort " + prothint.workDir + "/result_prosplign.gff)"
        diffResult += subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)

        if (diffResult == 0):
            shutil.rmtree(prothint.workDir + "/ProSplign")

    def testProSplignNonPBS(self):
        prothint.workDir = testDir + "/test_ProSplign/no_pbs"
        prothint.proteins = prothint.workDir + "/proteins.fasta"
        prothint.runProSplign(pbs=False)

        command = "diff <(sort " + prothint.workDir + "/ProSplign/scored_introns.gff) \
            <(sort " + prothint.workDir + "/result_scored_introns.gff)"
        diffResult = subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)

        command = "diff <(sort " + prothint.workDir + "/ProSplign/prosplign.gff) \
            <(sort " + prothint.workDir + "/result_prosplign.gff)"
        diffResult += subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)

        if (diffResult == 0):
            shutil.rmtree(prothint.workDir + "/ProSplign")


if __name__ == '__main__':
    testDir = os.path.abspath(os.path.dirname(__file__))
    sys.path.append(testDir + "/../bin")

    import prothint
    prothint.threads = str(multiprocessing.cpu_count())
    prothint.binDir = testDir + "/../bin"

    unittest.main()
