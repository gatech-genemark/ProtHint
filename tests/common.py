#!/usr/bin/env python3
# Author: Tomas Bruna

import unittest
import subprocess


def compareFiles(file1, file2):
    command = "diff <(sort " + file1 + " ) \
            <(sort " + file2 + " )"
    return subprocess.call(command, shell=True, executable='/bin/bash')
