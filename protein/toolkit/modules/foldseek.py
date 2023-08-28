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
import pandas as pd

class Foldseek:

    def __init__(self, pdbFi):
        self.pdbFi = pdbFi 
        
    def run_cmd(self,foldseek_db,outFi):
        cmd = f'/dat1/nbt2/miniconda3/bin/foldseek easy-search {self.pdbFi} {foldseek_db} {outFi} tmp2 --format-output "query,target,fident,alnlen,mismatch,gapopen,qlen,qstart,qend,tstart,tend,evalue,bits,prob,alntmscore,qtmscore,ttmscore,u,t"'
        print(cmd)
        run = subprocess.run(cmd, shell=True)
       
        return run.returncode
    
    def parse(self, outFi):
        dt = []
        with open(outFi, 'r') as f:
            for li in f:
                query,target,fident,alnlen,mismatch,gapopen,qlen,qstart,qend,tstart,tend,evalue,bits,prob,alntmscore,qtmscore,ttmscore,u,t = li.strip('\n').split('\t')
                if float(ttmscore) > 0.3:
                    dt.append({
                        "query":query,
                        # "target":target.split('.')[0],
                        "target":target,
                        "qlen": int(qlen),
                        "qstart": int(qstart),
                        "qend": int(qend),
                        "ttmscore":  round(float(ttmscore),2),
                        "qtmscore":round(float(qtmscore),2),
                        "prob": round(float(prob),2),
                        "evalue":round(float(evalue),2),
                        "u":u,
                        "t":t,
                    })
        return dt


    def format_ECOD(self,outFi, outJsonFi, ecodInfoFi):
        data = self.parse(outFi)
        webDt = []
        with open(ecodInfoFi, 'rb') as fe:
            ecod = pickle.load(fe)
            for item in data:
                uid = item['target'].split('.')[0]
                uid2,ecod_domain_id,pdb,chain,pdb_range,t_name,f_name = ecod[uid]
                item['pdb'] =pdb
                item['chain'] =chain
                item['ecod_domain_id'] = ecod_domain_id
                # item['pdb_range'] =pdb_range
                item['t_name'] = t_name.strip('"')
                item['f_name'] = f_name.strip('"')
                webDt.append(item)
        with  open(outJsonFi, 'w') as fo:
            json.dump(json.dumps(webDt), fo)

    def format_SCOP(self,outFi, outJsonFi, InfoFi):
        data = self.parse(outFi)
        dt = defaultdict(list)
        with open(InfoFi) as f:
            next(f)
            for li in f:
                cell = li.strip().split("\t")
                dt[cell[0]] = cell # cell[0]: uid, domain_id
        webDt = []
        for item in data:
            uid,uniprot,pdbid,chain,pdb_range,*_ = item['target'].split('.')
            uid = item['target'].split('.')[0]
            cell = dt[uid]
            item['pdbid'] =pdbid
            item['chain'] =chain
            item['ecod_domain_id'] = uid
            # item['pdb_range'] =pdb_range
            item['t_name'] = cell[-2]
            item['f_name'] = cell[-1]
            webDt.append(item)
         
        with open(outJsonFi, 'w') as fo:
            json.dump(json.dumps(webDt), fo)

    def format_pdbDB(self,outFi, outJsonFi, InfoFi):
        data = self.parse(outFi)
        dt = defaultdict(list)
        with open(InfoFi, 'r') as f:
            next(f); next(f)
            for li in f:
                cell = li.strip('\n').split('\t')
                dt[cell[0].lower() ] = cell
            
        webDt = []
        for  item in data:
            # uid,uniprot,pdbid,chain,pdb_range,*_ = item['target'].split('.')
            uid = item['target'].split('.')[0]
            item['pdbid'] =uid.upper()
            item['desc'] = ''
            if uid in dt:
                cell = dt[uid]
                print(uid, cell)
                item['desc'] = cell[3]
            webDt.append(item)
         
        with open(outJsonFi, 'w') as fo:
            json.dump(json.dumps(webDt), fo)
    def format_AFDB(self,outFi, outJsonFi, InfoFi):
        data = self.parse(outFi)
        dt = defaultdict(list)
        with open(InfoFi) as f:
            next(f)
            for li in f:
                cell = li.strip().split("\t")
                dt[cell[0]] = cell # cell[0]: uid, domain_id
        webDt = []
        for item in data:
            uid = item['target'].split('-')[1]
            if uid in dt:
                cell = dt[uid]
                accessions, entry_name, data_class, protein_name, gene_name, organism, taxonomy_id, sequence_length,*_ =cell
                item['pdbid'] =uid
                item['desc'] = protein_name
                item['entry_name'] = entry_name
                item['organism'] = organism
            webDt.append(item)
         
        with open(outJsonFi, 'w') as fo:
            json.dump(json.dumps(webDt), fo)
            
def run_annotate( pdbFi,outDir):
    foldseek = Foldseek(pdbFi)

    ecod_out = os.path.join(outDir, 'ecod.txt')
    ecod_json = os.path.join(outDir, 'ecod.json')
    ecod_foldseekDB = '/dat1/dat/db/foldseek/ecod/F70/foldseek/foldseek_ECOD_F70'
    ecodInfoFi = "/dat1/dat/db/ECOD/F70/ecod.F70.pickle"
    foldseek.run_cmd(ecod_foldseekDB,ecod_out)
    foldseek.format_ECOD(ecod_out, ecod_json,ecodInfoFi)

    scop_out = os.path.join(outDir, 'scop.txt')
    scop_json = os.path.join(outDir, 'scop.json')
    scop_foldseekDB = '/dat1/dat/db/foldseek/scop2/foldseek/foldseek_scopDomain'
    scopInfoFi = "/dat1/nbt2/proj/21-prot/dat/Scope2/scop-cla-latest.tab.txt"
    foldseek.run_cmd(scop_foldseekDB,scop_out)
    foldseek.format_SCOP(scop_out, scop_json, scopInfoFi)
    
    #TODO add CATH
    # cath_out = os.path.join(outDir, 'cath.txt')
    # cath_json = os.path.join(outDir, 'cath.json')
    # cath_foldseekDB = '/dat1/dat/db/cath/F70/foldseek/foldseek_cath_F70'
    # cathInfoFi = "/dat1/dat/db/cath/F70/cath.F70.pickle"
    # foldseek.run_cmd(cath_foldseekDB,cath_out)
    # foldseek.format_cath(cath_out, cath_json,cathInfoFi)

    pdbDB_out = os.path.join(outDir, 'pdbDB.txt')
    pdbDB_json = os.path.join(outDir, 'pdbDB.json')
    pdbDB_foldseekDB = '/dat1/dat/db/foldseek/pdbDB/foldseek_PDB'
    pdbDBInfoFi = "/dat1/nbt2/proj/21-prot/dat/pdb/derived_data/index/entries.idx"
    foldseek.run_cmd(pdbDB_foldseekDB,pdbDB_out)
    foldseek.format_pdbDB(pdbDB_out, pdbDB_json,pdbDBInfoFi)

    AFDB_out = os.path.join(outDir, 'AFDB.txt')
    AFDB_json = os.path.join(outDir, 'AFDB.json')
    AFDB_foldseekDB = '/dat1/dat/db/foldseek/AFDB/foldseek_PDB_AFDB'
    AFDBInfoFi = "/dat1/dat/db/uniprot/subset/afdb.info"
    foldseek.run_cmd(AFDB_foldseekDB,AFDB_out)
    foldseek.format_AFDB(AFDB_out, AFDB_json,AFDBInfoFi)

class Main:
    def run_ecod(self, pdbFi,outDir):
        run_annotate( pdbFi, outDir)

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

    def pickle_scope(Fi="/dat1/nbt2/proj/21-prot/dat/Scope2/scop-cla-latest.tab.txt", outFi="scope-cla.pickle"):
        dt = defaultdict(list)
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip().split("\t")
                dt[cell[0]] = cell
        return dt
        # fo =  open(outFi, 'wb')
        # pickle.dump(dt, fo)
        # fo.close()
                


if __name__ == '__main__':
    print('')
    fire.Fire(Main)