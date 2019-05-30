#!/usr/bin/env python
# Author: Tomas Bruna

import unittest
import sys
import os
import subprocess


class TestSpalnSelection(unittest.TestCase):

    def testProcessOutput(self):
        prothint.workDir = testDir + "/test_SpalnSelection"
        prothint.filterSpalnPairs(10)
        os.chdir(prothint.workDir)

        command = "diff <(sort Spaln/pairs_filtered.out) \
                  <(sort result_pairs_filtered.out)"
        diffResult = subprocess.call(command, shell=True, executable='/bin/bash')
        self.assertEqual(diffResult, 0)

        os.remove("Spaln/pairs_filtered.out")


if __name__ == '__main__':
    testDir = os.path.abspath(os.path.dirname(__file__))
    sys.path.append(testDir + "/../bin")

    import prothint
    prothint.binDir = testDir + "/../bin"

    unittest.main()
