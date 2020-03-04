#!/usr/bin/env python3
# Author: Tomas Bruna

import unittest
import os
import subprocess
import shutil
import common


class TestProtHint(unittest.TestCase):

    def testIter(self):
        os.chdir(testDir + "/test_iter")

        subprocess.call("../../bin/prothint.py genome.fasta proteins.fasta --geneSeeds genemark.gtf \
                        --prevGeneSeeds genemarkPrev.gtf --prevSpalnGff prevSpaln.gff --workdir output", shell=True)

        self.assertEqual(common.compareFiles("output/prothint.gff", "test_output/prothint.gff"), 0)
        self.assertEqual(common.compareFiles("output/evidence.gff", "test_output/evidence.gff"), 0)
        self.assertEqual(common.compareFiles("output/diamond/diamond.out", "test_output/diamond/diamond.out"), 0)
        self.assertEqual(common.compareFiles("output/uniqueSeeds.gtf", "test_output/uniqueSeeds.gtf"), 0)
        self.assertEqual(common.compareFiles("output/prevHints.gff", "test_output/prevHints.gff"), 0)


        shutil.rmtree("output")


if __name__ == '__main__':
    testDir = os.path.abspath(os.path.dirname(__file__))
    unittest.main()
