#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : structure_prediction.py
# @Date            : 2023/07/13 19:04:02
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip,subprocess


class RoseTTAFold2NA:

    def __init__(self,work_dir, ):
        self.work_dir = os.path.abspath(work_dir)
        self.outdir = os.path.join(self.work_dir, 'out')
        os.system('mkdir -p ' + self.outdir)
        self.init_cmd = 'ssh -p3389 nong@172.22.148.150 "source ~/.zshrc;cd {self.work_dir}; conda run -n RF2NA /dat1/apps/RoseTTAFoldNA/RoseTTAFold2NA/run_RF2NA.sh '
    def run(self,protein_fasta, nucleic_type, nucleic_fasta):
        cmd =self.init_cmd + f'{self.outdir} {protein_fasta} {nucleic_type}:{nucleic_fasta}"'
        subprocess.run(cmd, shell=True)
        self.clear()

    def clear(self,):
        cmd = 'cd {self.outdir};rm -rf db  hhblits log  RNA.afa  rna_binding_protein.atab  rna_binding_protein.hhr  rna_binding_protein.msa0.a3m  rna_binding_protein.RNA.a3m  RNA.fa  RNA.unfilter.afa  RNA.wquery.unfilt.afa  trim.db'
        subprocess.run(cmd, shell=True)

class Main:
    def RF2NA(self, work_dir,protein_fasta, nucleic_type, nucleic_fasta):
        rf2na = RoseTTAFold2NA(work_dir)
        rf2na.run(protein_fasta, nucleic_type, nucleic_fasta)

if __name__ == '__main__':
    fire.Fire(Main)