
import fire, os,subprocess
import pandas as pd

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

class Mymetal:

    def __init__(self,fasta):
        self.fasta = fasta
    
    def run(self):
        self.metal = mbp.predict(self.fasta)

    def save(self):
        iof.save_out_csv(self.metal, 'out_filename.csv')

def pdb2fasta(pdbFile, FaFi):
    if not os.path.exists(FaFi):
        cmd2 = f'pdb2fasta {pdbFile} > {FaFi}'
        subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)


def lMetalSite2protvista(score_, site_,label):
    fragments = []
    scores = score_.split(',')
    sites  = site_.split(',')

    for index, char in enumerate(sites):
        if char == "1":
            fragment_ = {
                "start": index + 1,
                "end": index + 1, 
                "tooltipContent": f"{label}, Probability: {scores[index]}", 
                "color":  map_color(float(scores[index]),0,1)
            }
            fragments.append(fragment_)
    
    data_subtrack = {
                    "accession": 'LMetalSite' ,
                    "type": 'Binding Site:',
                    "label": f"{label}" ,
                    "shape": 'circle',
                    "locations": [{"fragments":fragments}]
                    }
    
    return data_subtrack


class LMetalSite:
    def __init__(self, pdbFile, root_dir):
        self.pdbFile = os.path.abspath( pdbFile)
        self.root_dir = os.path.abspath(root_dir)
        self.FaFi = os.path.join(self.root_dir, 'seqs.fasta')
        self.outDir = os.path.join(self.root_dir, self.__class__.__name__ )
        self.outFi =  os.path.join(self.root_dir, self.__class__.__name__, 'seqs_predictions.csv')
        subprocess.run('mkdir -p ' + self.outDir, stdout=subprocess.PIPE, shell=True)
        self.res = []
    def run(self):
        # pdb2fasta(self.pdbFile, self.FaFi)
        cmd = f"conda run -n LMetalSite python /apps_dk/LMetalSite/script/LMetalSite_predict.py --fasta {self.FaFi} --outpath {self.outDir}"
       
        print(cmd)
        subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        self.parse()
    def parse(self):
        dt = pd.read_csv(self.outFi,dtype=str)
        print(dt.columns)
        ZN = lMetalSite2protvista(dt.ZN_prob[0], dt.ZN_pred[0], 'Zn')
        CA = lMetalSite2protvista(dt.CA_prob[0], dt.CA_pred[0], 'Ca')
        MG = lMetalSite2protvista(dt.MG_prob[0], dt.MG_pred[0], 'Mg')
        MN = lMetalSite2protvista(dt.MN_prob[0], dt.MN_pred[0], 'Mn')
        data_subtrack = []
        if len(ZN['locations'][0]['fragments']) >0:
            print("ZN:", ZN)
            data_subtrack.append(ZN)
        if len(CA['locations'][0]['fragments']) >0:
            data_subtrack.append(CA) 
        if len(MG['locations'][0]['fragments']) >0:
            data_subtrack.append(MG)
        if len(MN['locations'][0]['fragments']) >0:
            data_subtrack.append(MN)

        track = {
           "label": f'Metal Ion',
           "labelType": 'text',
           "data":data_subtrack
        }
        self.res =  track

def run_LMetalSite(pdbFile, root_dir):
    lms = LMetalSite(pdbFile, root_dir)
    lms.run()
    print(lms.res)
    return lms.res
class Main:
    def run(self,fasta):
        from mymetal import mbp,iof
        metal = Mymetal(fasta)
        metal.run()
        metal.save
    def run_lms(self,pdbFile,root_dir):
        run_LMetalSite(pdbFile,root_dir)

if __name__ == "__main__":
    fire.Fire(Main)