#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : foldseek.py
# @Date            : 2023/07/15 18:49:17
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip,subprocess,json
from collections import defaultdict
import pickle

class Foldseek:

    def __init__(self, pdbFi, foldseek_db,outFi):
        self.pdbFi = pdbFi 
        self.foldseek_db = foldseek_db
        self.outFi = outFi
        self.dt = []
    def run_cmd(self):
        cmd = f'foldseek easy-search {self.pdbFi} {self.foldseek_db} {self.outFi} tmp2 --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,alntmscore,qtmscore,ttmscore"'
        run = subprocess.run(cmd, shell=True)
        self.parse()
        return run.returncode
    def parse(self,):
        with open(self.outFi, 'r') as f:
            for li in f:
                query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,alntmscore,qtmscore,ttmscore = li.strip('\n').split('\t')
                self.dt.append({
                    "query":query,
                    "target":target.split('.')[0],
                })

    def format_ECOD(self,ecodFi,outJsonFi):
        webDt = []
        with open(ecodFi, 'rb') as fe:
            ecod = pickle.load(fe)
            for item in self.dt:
                uid = item['target']
                uid2,ecod_domain_id,pdb,chain,pdb_range,t_name,f_name = ecod[uid]
                item['pdb'] =pdb
                item['chain'] =chain
                item['pdb_range'] =pdb_range
                item['t_name'] = t_name
                item['f_name'] = f_name
                webDt.append(item)
        with  open(outJsonFi, 'w') as fo:
            json.dump(json.dumps(webDt), fo)



class Main:
    def run_ecod(self, pdbFi,foldseek_db,outFi):
        foldseek_ecod = Foldseek(pdbFi,foldseek_db,outFi)
        foldseek_ecod.run_cmd()
        foldseek_ecod.format_ECOD('/dat1/dat/db/ECOD/F70/ecod.F70.pickle','ecod.json')
    def pickle_ECOD(self,ecodFi="/dat1/dat/db/ECOD/F70/ecod.latest.F70.domains.txt",outFi="ecod.F70.pickle"):
        dt = defaultdict(list)
        with open(ecodFi, 'r') as f:
            for line in f:
                if re.match('#',line):continue
                cell = line.strip('\n').split('\t')
                uid,ecod_domain_id = cell[0:2]
                pdb,chain,pdb_range=cell[4:7]
                t_name,f_name = cell[12:14]
                dt[uid] = [uid,ecod_domain_id,pdb,chain,pdb_range,t_name,f_name]
        fo =  open(outFi, 'wb')
        pickle.dump(dt, fo)
        fo.close()


                


if __name__ == '__main__':
    print('')
    fire.Fire(Main)