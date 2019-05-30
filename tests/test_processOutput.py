#!/usr/bin/env python
# Author: Tomas Bruna

import unittest
import sys
import os
import filecmp


class TestProcessOutput(unittest.TestCase):

    def testProcessOutput(self):
        prothint.workDir = testDir + "/test_processOutput"
        prothint.processOutput(4, 0.3, 90)
        os.chdir(prothint.workDir)

        self.assertEqual(filecmp.cmp("introns.gff",
                                     "test_introns.gff"), True)

        self.assertEqual(filecmp.cmp("ProSplign/prosplign_combined.gff",
                                     "test_prosplign_combined.gff"), True)

        self.assertEqual(filecmp.cmp("hints.gff",
                                     "test_hints.gff"), True)

        self.assertEqual(filecmp.cmp("evidence.gff",
                                     "test_evidence.gff"), True)

        os.remove("introns.gff")
        os.remove("ProSplign/prosplign_combined.gff")
        os.remove("ProSplign/starts_stops.gff")
        os.remove("hints.gff")
        os.remove("evidence.gff")


if __name__ == '__main__':
    testDir = os.path.abspath(os.path.dirname(__file__))
    sys.path.append(testDir + "/../bin")

    import prothint
    prothint.binDir = testDir + "/../bin"

    unittest.main()
