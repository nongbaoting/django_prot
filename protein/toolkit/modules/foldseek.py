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

    def __init__(self, pdbFi):
        self.pdbFi = pdbFi 
        
    def run_cmd(self,foldseek_db,outFi):
        cmd = f'foldseek easy-search {self.pdbFi} {foldseek_db} {outFi} tmp2 --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,alntmscore,qtmscore,ttmscore"'
        run = subprocess.run(cmd, shell=True)
        return run.returncode
    
    def parse(self, outFi):
        dt = []
        with open(outFi, 'r') as f:
            for li in f:
                query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,alntmscore,qtmscore,ttmscore = li.strip('\n').split('\t')
                if float(ttmscore) > 0.3:
                    dt.append({
                        "query":query,
                        # "target":target.split('.')[0],
                        "target":target.split('.')[0],
                        "qstart": int(qstart),
                        "qend": int(qend),
                        "ttmscore":  round(float(ttmscore),2),
                        "qtmscore":round(float(qtmscore),2),
                        "prob": float(prob),
                        "evalue":float(evalue),
                    })
        return dt


    def format_ECOD(self,outFi, outJsonFi, ecodInfoFi):
        data = self.parse(outFi)
        webDt = []
        with open(ecodInfoFi, 'rb') as fe:
            ecod = pickle.load(fe)
            for item in data:
                uid = item['target']
                uid2,ecod_domain_id,pdb,chain,pdb_range,t_name,f_name = ecod[uid]
                item['pdb'] =pdb
                item['chain'] =chain
                item['ecod_domain_id'] = ecod_domain_id
                item['pdb_range'] =pdb_range
                item['t_name'] = t_name.strip('"')
                item['f_name'] = f_name.strip('"')
                webDt.append(item)
        with  open(outJsonFi, 'w') as fo:
            json.dump(json.dumps(webDt), fo)



def run_annotate( pdbFi,outDir):
    foldseek = Foldseek(pdbFi)
    ecod_out = os.path.join(outDir, 'ecod.txt')
    ecod_json = os.path.join(outDir, 'ecod.json')
    ecod_foldseekDB = '/dat1/dat/db/ECOD/F70/foldseek/foldseek_ECOD_F70'
    ecodInfoFi = "/dat1/dat/db/ECOD/F70/ecod.F70.pickle"
    #foldseek.run_cmd(ecod_foldseekDB,ecod_out)
    foldseek.format_ECOD(ecod_out, ecod_json,ecodInfoFi)

class Main:
    def run_ecod(self, pdbFi,outDir):
        run_annotate( pdbFi,outDir)

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