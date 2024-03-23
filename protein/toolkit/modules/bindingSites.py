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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
def map_color(value, min_value, max_value):
    normalized_value = (value - min_value) / (max_value - min_value)
    cmap = plt.cm.get_cmap('RdYlGn')  # 选择色谱，这里使用红黄绿色谱
    rgba = cmap(normalized_value)
    red = int(rgba[0] * 255)
    green = int(rgba[1] * 255)
    blue = int(rgba[2] * 255)
    color = (red, green, blue)
    hex_color = '#%02x%02x%02x' % color
    return hex_color

# min_value = float(input("请输入最小值："))
# max_value = float(input("请输入最大值："))
# number = float(input("请输入一个数字："))

# color = map_color(number, min_value, max_value)
# print("对应的颜色为：", color)


def color_pick(c):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcol
    import matplotlib.cm as cm
    start_time = 0
    end_time = 1

    # Generate some dummy data.
    tim = range(start_time,end_time)

    # Make a user-defined colormap.
    cm_br = mcol.LinearSegmentedColormap.from_list("MyCmapName",["b","r"])

    # Make a normalizer that will map the time values from
    # [start_time,end_time+1] -> [0,1].
    cnorm = mcol.Normalize(vmin=min(tim),vmax=max(tim))

    # Turn these into an object that can be used to map time values to colors and
    # can be passed to plt.colorbar().
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm_br)
    cpick.set_array([])
    a = cpick.to_rgba(c)

    return mcol.to_hex(a)

class ScanNet:

    def __init__(self,pdbFile,root_dir):
        self.init_cmd = 'conda run -n py_scannet python /apps_dk/ScanNet/predict_bindingsites.py '
        self.root_dir = os.path.abspath(root_dir)
        self.pdbFile = os.path.abspath( pdbFile)
        self.name = os.path.basename(pdbFile).split('.')[0]
    
    def run(self):
        cmd = self.init_cmd + f'{self.pdbFile}  --noMSA --predictions_folder {self.root_dir}/scannet_out'
        print(cmd)
        subprocess.run(cmd, shell=True)

    def parse2protvista(self):
        outdir = f'{self.root_dir}/scannet_out/{self.name}_single_ScanNet_interface_noMSA/'
        outCsv = outdir + f'predictions_{self.name}.csv'
        outPdb = outdir + f'annotated_{self.name}.pdb'
        track = {
           "label": 'PPI Sites',
           "labelType": 'text',
           "data":[]

        }
        fragments = []

        with open(outCsv,'r') as f:
            next(f)
            for index, line in enumerate(f):
                Model,Chain,Residue_Index,Sequence,prob = line.strip().split(',')
                if float(prob) >0.5:
                    fragment_ = {
                        "start": index + 1,
                        "end": index + 1, 
                        "tooltipContent": f"Position: {Residue_Index}, Seq: {Sequence}, Score: {prob}", 
                        "color": map_color(float(prob),0,1)
                    }
                    fragments.append(fragment_)


        fragments_dt = {"fragments":fragments}
        data_subtrack = {
                        "accession": 'ScanNet' ,
                        "type": 'PPI',
                        "label": "ScanNet" ,
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
    def run_scannet(self,pdbFile,root_dir):
        import json
        track = run_ppi(pdbFile, root_dir)
        with open('track.json', 'w') as fp:
            json.dump(json.dumps(track), fp)


class CLAPE:
    def __init__(self,pdbFile, root_dir):
        self.init_cmd = 'ssh -p3389 nong@172.22.148.150 "source ~/.zshrc;'
        self.root_dir = os.path.abspath(root_dir)
        self.pdbFile = os.path.abspath( pdbFile)
        self.name = os.path.basename(pdbFile).split('.')[0]

    def run(self):
        cmd = "conda run -n CLAPE python /dat1/apps2/CLAPE/clape.py --input example.fa --output CLAPE.txt"

if __name__ == '__main__':
    print('')
    fire.Fire(Main)