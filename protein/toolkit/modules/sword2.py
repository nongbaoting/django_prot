#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : sword2.py
# @Date            : 2023/04/20 21:32:28
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip,json
import subprocess
from Bio import SeqIO

colorSet3D = [
"#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", 
"#666666", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", 
"#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
 "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
]

def scanAndFind_pattern(mydir, mypattern):
    wantFiles = []
    for entry in os.scandir(mydir):
        if entry.is_file() and mypattern.search(entry.name):
            wantFiles.append(entry)
        elif entry.is_dir():
            wantFiles.extend(scanAndFind_pattern(entry.path, mypattern))
    return wantFiles

def domain2track(domain_index,domains):
    domain_track = {
        "trackId": f"track_{domain_index}",
        "rowTitle": f"Architecture_{domain_index}",
        "titleFlagColor":colorSet3D[domain_index],
        "displayType": 'composite',
        "trackHeight": 25,
        "displayConfig": [

        ]
    }
    re_cc = re.compile(r';')
    re_seg = re.compile(r'-')
    for i,domain in enumerate(domains):
        units = domain.split(';')
        locs = [i.split("-") for i in units]
        locs_flat = [int(item.strip()) for sublocs in locs for item in sublocs]
        domain_start,domain_end = min(locs_flat), max(locs_flat)
        domain_sub = re_cc.sub('/', domain)
        domain_sub =re_seg.sub('~', domain_sub)
        displayConfig = {
            "displayType": 'block',
            "displayColor": colorSet3D[i],
            "displayData": [
                {
                "begin": domain_start,
                "end": domain_end,
                "gaps": [
                
                ],
                "displayId": f'domain_{i}',
                "sourceId": colorSet3D[i],
                "featureId": domain_sub,
                },
            ],
        }
        if len(units)>1:
            unit_loc = [[int(ii.strip()) for ii in i.split('-')] for i in units]
            print(unit_loc)
            unit_loc = sorted(unit_loc, key=lambda k:k[0])
            print(unit_loc)
            gaps = []
            for gap_i in range(len(unit_loc) -1 ):
                gap_start = unit_loc[gap_i][1] + 1
                gap_end = unit_loc[gap_i+1][0] -1
                gap ={
                    "begin": gap_start, "end": gap_end 
                }
                gaps.append(gap)
            displayConfig["displayData"][0]["gaps"] =gaps
        domain_track["displayConfig"].append(displayConfig)
    return domain_track

def domain2protvistaFragments(domain_index,domains):
    re_cc = re.compile(r';')
    re_seg = re.compile(r'-')
    fragments = []

    for i,domain in enumerate(domains):
        units = domain.split(';')
        for domain_units in units:
            locs = [int(ii.strip()) for ii in domain_units.split('-')]
            domain_start,domain_end = locs
            fragment_ = {
                "start": domain_start,
                "end": domain_end, 
                "tooltipContent": f"Type: domain {i}<br>Range: {domain_start} - {domain_end}", 
                "color": colorSet3D[i], 
                }

            fragments.append(fragment_)
    return {"fragments":fragments}


class Sword2:
    def __init__(self,):
        self.suffix = 'sword.txt'
        self.suffix_reg=  re.compile(f"{self.suffix}$")
        self.sep_reg = re.compile(f"/")
        self.fasta = {}
    
    def scandir(self,mydir,outFi,outFiProtvista):
        entries = scanAndFind_pattern(mydir, self.suffix_reg)
        dt = {}
        dt_pt = {}
        for entry in entries:
            entry_id = self.sep_reg.split(entry.path)[-2].split('_')[0]
            seq_id = entry_id 
            data = self.parse(entry.path, seq_id)
            data_pt = self.parse2protvista(entry.path, seq_id)
            dt[entry_id] = data
            dt_pt[entry_id] = data_pt
        fo = open(outFi,'w')
        json.dump(json.dumps(dt), fo)
        fo.close()
        with open(outFiProtvista, 'w') as fp:
            json.dump(json.dumps(dt_pt), fp)
        return dt_pt
    def parse(self, Fi, seq_id):
        seq  = str(self.fasta[seq_id].seq)
        rowConfigData = [
            {
            "trackId": 'sequenceTrack',
            "trackHeight": 20,
            "trackColor": '#FEFEFE',
            "displayType": 'sequence',
            "nonEmptyDisplay": True,
            "rowTitle": 'SEQUENCE',
            "trackData": [
                {
                "begin": 1,
                "value":seq,
                },
                ],
                }
            ]
        re_numbers = re.compile(r'^\d')
        with open(Fi, "r") as f:
            num = 0
            for line in f:
                if re_numbers.match(line):
                    print(line)
                    domain_num,min_len,boundary,ave_k,quality,_ = line.strip('\n').split('|')
                    domains = boundary.strip().split(' ')
                    track = domain2track(num,domains)
                    num +=1
                    rowConfigData.append(track)
                    data= {
                        "sequence": seq,
                        "rowConfigData": rowConfigData
                    }
        return data
    
    def parse2protvista(self, Fi, seq_id):
        seq  = str(self.fasta[seq_id].seq)
        track = {
           "label": 'Sword2 Domains',
           "labelType": 'text',
           "data":[]

        }
        re_numbers = re.compile(r'^\d')
        with open(Fi, "r") as f:
            num = 0
            for line in f:
                if re_numbers.match(line):
                    print(line)
                    domain_num,min_len,boundary,ave_k,quality,_ = line.strip('\n').split('|')
                    domains = boundary.strip().split(' ')
                    fragments = domain2protvistaFragments(num,domains)
                    data_subtrack = {
                        "accession": 'sword' + str(num),
                        "type": 'Sword2 Domain',
                        "label": "Archi. " + str(num),
                        "labelTooltip":"Sword2 Archi " + str(num),
                        "locations": [fragments]
                    }
                    num +=1
                    track['data'].append(data_subtrack)
                    
        return track


    def parseFasta(self,FaFi):
        self.fasta = SeqIO.to_dict(SeqIO.parse(FaFi, 'fasta'))

    def run(self,pdbfile,chain,work_dir):
        
        cmd = f'cd {work_dir};conda run -n sword2 /apps_dk/SWORD2/SWORD2.py -i  {pdbfile} -o result -c {chain} -x 8'
        
        subprocess.run(cmd, shell=True)

def run_sword2(pdbfile,chain, work_dir):
    sword2 = Sword2()
    fasta_seq = os.path.join(work_dir,'seqs.fasta')
    pdbfile_ln = os.path.join(work_dir,'upload.pdb')
    cmd = f'cp  {pdbfile} {pdbfile_ln}'
    subprocess.run(cmd,shell=True)

    sword2.run(pdbfile_ln,chain,work_dir)
    sword2.parseFasta(fasta_seq)
    outFi = os.path.join(work_dir,'sword2.track.json')
    outFiProtvista = os.path.join(work_dir,'sword2.protvistaTrack.json')
    sword2_track = sword2.scandir(work_dir,outFi,outFiProtvista)
    seq_id = list(sword2.fasta.keys())[0]
    return  sword2_track,str(sword2.fasta[seq_id].seq)

class Main:
    def run(self,pdbfile,chain,work_dir):
        run_sword2(pdbfile, chain, work_dir)

if __name__ == '__main__':
    print('')
    fire.Fire(Main)