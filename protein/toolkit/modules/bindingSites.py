#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : bindingSites.py
# @Date            : 2023/07/10 19:26:03
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip,subprocess
# from functions import *
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


class ScanNet:

    def __init__(self,pdbFile,root_dir):
        self.init_cmd = 'ssh -p3389 nong@172.22.148.150 "source ~/.zshrc;conda run -n py_scannet python /dat1/apps/ScanNet/predict_bindingsites.py '
        self.root_dir = os.path.abspath(root_dir)
        self.pdbFile = os.path.abspath( pdbFile)
        self.name = os.path.basename(pdbFile).split('.')[0]
    
    def run(self):
        cmd = self.init_cmd + f'{self.pdbFile}  --noMSA --predictions_folder {self.root_dir}/scannet_out"'
        print(cmd)
        subprocess.run(cmd, shell=True)


    def parse2protvista(self):
        outdir = f'{self.root_dir}/scannet_out/{self.name}_single_ScanNet_interface_noMSA/'
        outCsv = outdir + f'predictions_{self.name}.csv'
        outPdb = outdir + f'annotated_{self.name}.pdb'
        track = {
           "label": 'Protein Binding Sites',
           "labelType": 'text',
           "data":[]

        }
        fragments = []

        with open(outCsv,'r') as f:
            next(f)
            for line in f:
                Model,Chain,Residue_Index,Sequence,prob = line.strip().split(',')
                if float(prob) >0.5:
                    fragment_ = {
                        "start": int(Residue_Index),
                        "end": int(Residue_Index), 
                        "tooltipContent": f"AA: {Sequence},Score: {prob}", 
                        "color": colorSet3D[0], 
                    }
                    fragments.append(fragment_)


        fragments_dt = {"fragments":fragments}
        data_subtrack = {
                        "accession": 'ScanNet' ,
                        "type": 'ScanNet',
                        "label": "ScanNet " ,
                         "shape": 'circle',
                        "labelTooltip":"ScanNet PPI sites" ,
                        "locations": [fragments_dt]
                    }
        track = {
           "label": 'Protein Binding Sites',
           "labelType": 'text',
           "data":[data_subtrack]

        }
        return track

def run_ppi(pdbFile,root_dir):
    scannet = ScanNet(pdbFile,root_dir)
    scannet.run()
    return scannet.parse2protvista()
class Main:
    def run(self,pdbFile,root_dir):
        import json
        track = run_ppi(pdbFile,root_dir)
        with open('track.json', 'w') as fp:
            json.dump(json.dumps(track), fp)

if __name__ == '__main__':
    print('')
    fire.Fire(Main)