#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : unidoc.py
# @Date            : 2023/04/20 10:27:09
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip,json,subprocess
from Bio import SeqIO
colorSet3D = [
    "#27A3B4",
    # "#C08423",
    "#D83E7C",
    "#1B9E77",
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#1B9E77',
    '#D95F02',
    '#7570B3',
    '#E7298A',
    '#66A61E',
    '#E6AB02',
    '#A6761D',
]

def scanAndFind_pattern(mydir, mypattern):
    wantFiles = []
    for entry in os.scandir(mydir):
        if entry.is_file() and mypattern.search(entry.name):
            wantFiles.append(entry)
        elif entry.is_dir():
            wantFiles.extend(scanAndFind_pattern(entry.path, mypattern))
    return wantFiles

def domain2track(domains):
    domain_track = {
        "trackId": "Domain_track",
        "rowTitle": "Unidoc DOMAIN",
        "titleFlagColor":"#27A3B4",
        "displayType": 'composite',
        "trackHeight": 25,
        "displayConfig": [
        ]
    }
    for i,domain in enumerate(domains):
        units = domain.split('/')
        locs = [i.split("~") for i in units]
        locs_flat = [int(item) for sublocs in locs for item in sublocs]
        domain_start,domain_end = min(locs_flat), max(locs_flat)
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
                "featureId": domain,
                },
            ],
        }
        if len(units)>1:
            unit_loc = [[int(ii) for ii in i.split('~')] for i in units]
            unit_loc = sorted(unit_loc, key=lambda k:k[0])
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
    # print(domain_track)
    return domain_track



class Unidoc:

    def __init__(self,pdbfile,chain,work_dir):
        self.pdbfile = pdbfile
        self.chain = chain
        self.work_dir = work_dir
        self.fasta = {}
        self.unidoc_out = os.path.join(self.work_dir,'res.unidoc')
        self.FaFi = os.path.join(self.work_dir,'seqs.fasta')

    def parse_unidoc(self):
        dt = {}
        sep_reg = re.compile(f".pdb|.cif")
        fileName = os.path.basename(self.pdbfile)
        entry_id = sep_reg.split(fileName)[0]
        seq_id = fileName.split('.')[0] + f':{self.chain}'
        seq  = str(self.fasta[seq_id].seq)
        domains = open(self.unidoc_out).read().strip().split(',')
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
        track = domain2track(domains)
        # print(track)
        rowConfigData.append(track)
        data= {
        "sequence": seq,
        "rowConfigData": rowConfigData
        }
        dt[entry_id] = data
        
        fo = open(f'{self.work_dir}/unidoc.track.json','w')
        json.dump(json.dumps(dt), fo)
    def parse2protvista(self,):
        dt = {}
        sep_reg = re.compile(f".pdb|.cif")
        fileName = os.path.basename(self.pdbfile)
        entry_id = sep_reg.split(fileName)[0]
        seq_id = fileName.split('.')[0] + f':{self.chain}'
        seq  = str(self.fasta[seq_id].seq)
        domains = open(self.unidoc_out).read().strip().split(',')
        
        # fragments is domain loctions
        fragments = []
        d_num = 0
        for i,domain in enumerate(domains):
            units = domain.split('/')
            for j, domain_units in enumerate(units):
                locs = [int(ii.strip()) for ii in domain_units.split('~')]
                domain_start,domain_end = locs
                fragment_ = {
                    "start": domain_start,
                    "end": domain_end, 
                    "tooltipContent": f"Type: domain {i}:{j}<br>Range: {domain_start} - {domain_end}", 
                    "color": colorSet3D[d_num], 
                    }
                d_num +=1
                fragments.append(fragment_)

            
        data_subtrack = {
                        "accession": 'Unidoc' ,
                        "type": 'Unidoc Domain',
                        "label": "Unidoc",
                        "labelTooltip":"Unidoc Domain",
                        "locations": [{"fragments":fragments}
                            
                        ]
                    }
        track = {
           "label": 'Unidoc Domains',
           "labelType": 'text',
           "data":[data_subtrack]

        }
        return track

    def parseFasta(self):
        self.fasta = SeqIO.to_dict(SeqIO.parse(self.FaFi, 'fasta'))

    def run(self):
        cmd1 = f'python /dat1/apps/structure/UniDoc/Run_UniDoc_from_scratch_structure.py -i {self.pdbfile}  -c {self.chain} -o {self.unidoc_out}'
        cmd2 = f'pdb2fasta {self.pdbfile} > {self.FaFi}'
        print(cmd1)
        subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)
        subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)
         
def run_unidoc(pdbfile, chain, work_dir):
    unidoc = Unidoc(pdbfile,chain,work_dir)
    unidoc.run()
    unidoc.parseFasta()
    unidoc.parse_unidoc()
    return unidoc.parse2protvista()

class Main:
    def run(self,pdbfile, chain, work_dir):
        run_unidoc(pdbfile, chain, work_dir)

if __name__ == '__main__':
    fire.Fire(Main)