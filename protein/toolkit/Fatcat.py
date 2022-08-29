
import re
import fire
import os
from collections import defaultdict
from crispr.utils.myFunctions import scanAndFind_pattern, calAlignIdent, run_cmd
from Bio import PDB
from protein.toolkit import myFunctions


class Fatcat:

    def __init__(self, path_1, path_2, tempDir, outFi):
        self.path_1 = path_1
        self.path_2 = path_2
        self.tempDir = tempDir
        self.matrix = []
        self.logFi = ''
        self.outFi = outFi

    def align(self):
        os.system("mkdir -p " + self.tempDir)
        dirname_1 = os.path.dirname(self.path_1)
        dirname_2 = os.path.dirname(self.path_2)
        basename_1 = os.path.basename(self.path_1)
        basename_2 = os.path.basename(self.path_2)

        cmd = f'cd {self.tempDir}; FATCAT  -i1 {dirname_1} -i2 {dirname_2} -p1 {basename_1} -p2 {basename_2} -q -t; cp tmp.opt.twist.pdb { self.outFi }'

        run = run_cmd(cmd)
        log = run.stdout.decode('utf-8')
        self.content = {
            "log": log
        }

    def save_pickle(self,):
        item_pickle = self.outFi + '.pickle'
        myFunctions.pickle_dumpObj2file(self.content, item_pickle)


def pair_align(path_1, path_2, tempDir, outFi_name):
    outFi = os.path.join(tempDir, outFi_name)
    fatcat = Fatcat(path_1, path_2, tempDir, outFi)
    fatcat.align()
    fatcat.save_pickle()
    return fatcat


if __name__ == "__main__":
    fire.Fire(pair_align)
