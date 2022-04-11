#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : TMalign.py
# @Date            : 2022/03/14 22:35:35
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os
from os import path
import re
import sys
import fire
import subprocess
from multiprocessing import Pool
from collections import defaultdict
from protein.toolkit import myFunctions

def parse_TMalign(outLogs):
    re_pdb = re.compile(r'.pdb$|.pdb.gz$|.cif$|.cif.gz$')
    re_TMalign_sp = re.compile(r"TM-align \(Version")
    re_TMalign_nameLen = re.compile(
        r"Name of Chain_1: (?P<chain_1>.*)\(to be superimposed onto Chain_2\)\nName of Chain_2:\s+(?P<chain_2>.*)\nLength of Chain_1:\s+(?P<chain_1_len>\d+) residues\nLength of Chain_2:\s+(?P<chain_2_len>\d+)\s+residues\n", re.M)
    re_TMalign_misc = re.compile(
        r"Aligned length=\s+(?P<align_len>\d+), RMSD=\s+(?P<RMSD>\d+(\.\d+)?), Seq_ID=n_identical/n_aligned=\s+(?P<Seq_ID>\d+(\.\d+)?)", re.M)
    re_TMalign_score = re.compile(
        r"TM-score=\s+(?P<tmscore_1>\d(\.\d+)?) \(if normalized by length of Chain_1, i.e., LN=\d+, d0=(?P<d01>\d+(\.\d+)?)\)\nTM-score= (?P<tmscore_2>\d(\.\d+)?) \(if normalized by length of Chain_2, i.e., LN=\d+, d0=(?P<d02>\d+(\.\d+)?)\)", re.M)
    re_TMalign_align = re.compile(
        r'\(":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues\)\n(?P<seq_1>.*)\n(?P<pairwise>.*)\n(?P<seq_2>.*)\n\nTotal CPU time')
    # outLogs = open(logFile, 'r').read().strip()
    lt = []
    for report in re_TMalign_sp.split(outLogs)[1:]:
        m_nameLen = re_TMalign_nameLen.search(report)
        m_misc = re_TMalign_misc.search(report)
        m_score = re_TMalign_score.search(report)
        m_align = re_TMalign_align.search(report)
        # print(m_nameLen.group())
        chain_1 = path.basename(m_nameLen.group('chain_1').strip())
        chain_2 = path.basename(m_nameLen.group('chain_2').strip())
        # chain_1 = re_pdb.sub('', chain_1)
        # chain_2 = re_pdb.sub('', chain_2)
        if not m_misc:
            print('m_misc')
            print(report)

        tmscore_1 = m_score.group('tmscore_1')
        tmscore_2 = m_score.group('tmscore_2')
        seq_1 = m_align.group('seq_1')
        seq_2 = m_align.group('seq_2')
        pairwise = m_align.group('pairwise')
        chain_1_len = m_nameLen.group('chain_1_len')
        chain_2_len = m_nameLen.group('chain_2_len')
        align_len = m_misc.group('align_len')
        cov_1 = round(float(align_len) / float(chain_1_len), 2)
        cov_2 = round(float(align_len) / float(chain_2_len), 2)
        lt_ = [chain_1, chain_2,
               m_nameLen.group('chain_1_len'), m_nameLen.group(
                   'chain_2_len'), m_misc.group('align_len'), str(cov_1), str(cov_2),
               m_misc.group('RMSD'), m_misc.group(
                   'Seq_ID'), tmscore_1, tmscore_2,
               m_score.group("d01"), m_score.group("d02"),
               m_align.group("seq_1"), m_align.group("pairwise"), m_align.group("seq_2")]
        lt_ = {
            "chain_1": chain_1, "chain_2":chain_2,
             "chain_1_len":  m_nameLen.group('chain_1_len'),  'chain_2_len': m_nameLen.group(
                   'chain_2_len'), 'align_len': m_misc.group('align_len'), 
                   "cov_1":str(cov_1), "cov_2": str(cov_2),
               'RMSD':m_misc.group('RMSD'), 'Seq_ID': m_misc.group(
                   'Seq_ID'), "tmscore_1": tmscore_1, "tmscore_2": tmscore_2,
               "d01": m_score.group("d01"), "d02":m_score.group("d02"),
               "seq_1": m_align.group("seq_1"), "pairwise": m_align.group("pairwise"), "seq_2": m_align.group("seq_2")
        }
        if float(tmscore_2) > 0.5:
            # print(lt_)
            # print("%s\n%s\n%s\n" % (seq_1,pairwise,seq_2))
            lt.append(lt_)
    return lt


def scanAndFind_pattern(mydir, mypattern):
    wantFiles = []
    for entry in os.scandir(mydir):
        if (entry.is_file() or entry.is_link()) and mypattern.search(entry.name):
            wantFiles.append(entry)
        elif entry.is_dir():
            wantFiles.extend(scanAndFind_pattern(entry.path, mypattern))
    return wantFiles

def run_cmd(cmd):
    run = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return run

def parse_uniprot(Fi='/dat1/nbt2/proj/21-prot/alphafold/uniprot_info/uniprot-proteome_UP000005640.tab'):
    dt = defaultdict(list)
    with open(Fi) as f:
        next(f)
        for li in f:
            cell = li.strip().split("\t")
            uniprot = cell[0]
            uniprot_name = cell[1]
            status = cell[2]
            protein_names = cell[3]
            gene_names = cell[4]
            organism = cell[5]
            length = cell[6]
            dt[uniprot] = [protein_names, gene_names]
    return dt


def parse_scope(Fi="/dat1/nbt2/proj/21-prot/dat/Scope2/scop-cla-latest.tab.txt"):
    dt = defaultdict(list)
    with open(Fi) as f:
        next(f)
        for li in f:
            cell = li.strip().split("\t")
            dt[cell[0]] = cell
    return dt

# scopeDomain_dir = "/dat1/nbt2/proj/21-prot/dat/pdb/scope_domain"
scopeDomain_dir ="/dat1/nbt2/proj/21-prot/dat/pdb/test"
class TMalgin:

    def __init__(self,params):
        self.dir_1 = params['dir_1']
        self.dir_2 = scopeDomain_dir
        self.result_dir = myFunctions.create_tmpDir(params['struc_cpm_dir'], params['uuid'])
        self.outFi = os.path.join(self.result_dir, params['outFi_name'])

    def run(self,):
        self.dir_vs_dir(self.dir_1, self.dir_2)
        self.parseLogFile()
        subprocess.run(f'cp  {self.dir_1}/* {self.result_dir}' , shell=True)
        return self.outFi

    def dir_vs_dir(self, dir1, dir2, core=16):
        pattern = re.compile('.pdb$|.pdb.gz$', re.IGNORECASE)
        pdbfiles_1 = scanAndFind_pattern(dir1, pattern)
        pdbfiles_2 = scanAndFind_pattern(dir2, pattern)
        # fo = open(outFi, 'w')
        res = []
        pool = Pool(core)
        for entry1 in pdbfiles_1:
            for entry2 in pdbfiles_2:
                cmd = f'TMalign {entry1.path} {entry2.path}'
                # cmd = f'SPalignNS -pair {entry1.path} {entry2.path}'
                res_ = pool.apply_async(run_cmd, (cmd,))
                # run = res_.get()
                # fo.write(run.stdout.decode('utf-8'))
                res.append(res_)
        pool.close()
        pool.join()
        self.outLogs = ''
        for res_ in res:
            run = res_.get()
            self.outLogs += run.stdout.decode('utf-8')
        #     fo.write(run.stdout.decode('utf-8'))
        # fo.close()

    def one_vs_all(self, input_pdb, db_dir, core=8):
        pattern = re.compile('.pdb$|.pdb.gz$', re.IGNORECASE)
        pdbfiles = scanAndFind_pattern(db_dir, pattern)

        res = []
        pool = Pool(core)
        for entry in pdbfiles:
            cmd = f'TMalign {input_pdb} {entry.path}'
            res_ = pool.apply_async(run_cmd, (cmd,))
        pool.close()
        pool.join()
        fo = open("res.txt", 'w')
        for res_ in res:
            run = res_.get()
            fo.write(run.stdout.decode('utf-8'))
        fo.close()

    def parseLogFile(self, ):
        tmalign = parse_TMalign(self.outLogs)
       
        heads = ['chain_1', 'chain_2', 'chain_1_len', 'chain_2_len', 'align_len', 'cov_1', 'cov_2',
                 'RMSD', 'Seq_ID', 'TMscore_1', 'TMscore_2', 'd0_1', 'd0_2',
                 'seq_1', 'pairwise', 'seq_2',
                 'chain_1_uniprot',
                 'chain_2_scopeDomain', 'chain_2_uniprot', 'chain_2_pdb', 'chain_2_chain'
                 ]
        # fo.write("\t".join(heads) + '\n')
        print("#format scope")
        arr = []
        for item in tmalign:
            pdb1 = item['chain_1']
            pdb2 = item['chain_2']
            # pdb1, pdb2 = arr[0:2]
            pdb1_prot = pdb1.split('.')[0]
            pdb2_domainID, pdb2_prot, pdb2_pdbID, chain = pdb2.split('.')[0:4]
            # arr.append(pdb1_prot)
            # arr.extend([pdb2_domainID, pdb2_prot, pdb, chain])
            # fo.write("\t".join(arr) + '\n')
            item['chain_1_uniprot'] = pdb1_prot

            item['chain_2_scopeDomain'] = pdb2_domainID
            item['chain_2_uniprot'] =  pdb2_prot
            item['chain_2_pdb'] = pdb2_pdbID
            item[ 'chain_2_chain'] = chain 
            arr.append(item)
        myFunctions.pickle_dump2file(arr, self.outFi)
      


if __name__ == "__main__":
    fire.Fire(TMalgin)
