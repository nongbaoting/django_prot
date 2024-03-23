
import os,subprocess
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

def musiteDeep2protvista(res):
    Phosphorylation,Glycosylation,Ubiquitination,SUMOylation,Acetylation,Methylation ,Pyrrolidone_carboxylic_acid,Palmitoylation,Hydroxylation= [],[],[],[],[],[],[],[],[]
    for res_ in res:
        pos, residue, ptm_type, score = res_
        fragment_ = {
            "start": pos,
            "end": pos, 
            "tooltipContent": f"{ptm_type}, Residue: {residue}, Score: {score}",
            "color":  map_color(float(score),0,1)
            }

        if  ptm_type == "Phosphotyrosine" or  ptm_type == "Phosphoserine_Phosphothreonine":
            Phosphorylation.append(ptm_type)
        elif ptm_type == "O-linked_glycosylation" or ptm_type == "N-linked_glycosylation":
            Glycosylation.append( fragment_ )
        elif ptm_type == "Methyllysine" or ptm_type == "Methylarginine":
            Methylation.append( fragment_ )
        elif ptm_type == "Ubiquitination":
            Ubiquitination.append( fragment_ )
        elif ptm_type == "SUMOylation":
            SUMOylation.append( fragment_ )
        elif ptm_type == "N6-acetyllysine":
            Acetylation.append( fragment_ )
        elif ptm_type == "Pyrrolidone_carboxylic_acid":
            Pyrrolidone_carboxylic_acid.append( fragment_ )
        elif ptm_type == "S-palmitoyl_cysteine":
            Palmitoylation.append( fragment_ )
        elif ptm_type == "Hydroxylysine" or ptm_type == "Hydroxyproline":
            Hydroxylation.append( fragment_ )

    data_subtrack_p = {
                "accession": 'PTM_me' ,
                "label": "Phosphorylation" ,
                "type": 'Binding Site:',
                "shape": 'circle',
                "locations": [{"fragments":Phosphorylation}]
                }
    data_subtrack_ub = {
                    "accession": 'PTM_ub' ,
                    "label": "Ubiquitination" ,
                     "type": 'Binding Site:',
                    "shape": 'circle',
                    "locations": [{"fragments":Ubiquitination}]
                    }
    data_subtrack_me = {
                    "accession": 'PTM_me' ,
                    "label": "Methylation" ,
                     "type": 'Binding Site:',
                    "shape": 'circle',
                    "locations": [{"fragments":Methylation}]
                    }
    data_subtrack_g = {
                    "accession": 'PTM_me' ,
                    "label": "Glycosylation" ,
                     "type": 'Binding Site:',
                    "shape": 'circle',
                    "locations": [{"fragments":Glycosylation}]
                    }
    data_subtrack_su = {
                    "accession": 'PTM_me' ,
                    "label": "SUMOylation" ,
                     "type": 'Binding Site:',
                    "shape": 'circle',
                    "locations": [{"fragments":SUMOylation}]
                    }
    data_subtrack_ac = {
                    "accession": 'PTM_me' ,
                    "label": "Acetylation" ,
                     "type": 'Binding Site:',
                    "shape": 'circle',
                    "locations": [{"fragments":Acetylation}]
                    }
    data_subtrack_pc = {
                    "accession": 'PTM_me' ,
                    "label": "Pyrrolidone carboxylic acid" ,
                     "type": 'Binding Site:',
                    "shape": 'circle',
                    "locations": [{"fragments":Pyrrolidone_carboxylic_acid}]
                    }
    data_subtrack_pa = {
                    "accession": 'PTM_me' ,
                    "label": "Palmitoylation" ,
                     "type": 'Binding Site:',
                    "shape": 'circle',
                    "locations": [{"fragments":Palmitoylation}]
                    }
    data_subtrack_h = {
                    "accession": 'PTM_me' ,
                    "label": "Hydroxylation" ,
                     "type": 'Binding Site:',
                    "shape": 'circle',
                    "locations": [{"fragments":Hydroxylation}]
                    }


    data_subtrack = []
    if len(Ubiquitination)  >0:
        data_subtrack.append(data_subtrack_ub)
    if len(Methylation)    >0:
        data_subtrack.append(data_subtrack_me)
    if len(Phosphorylation) >0:
        data_subtrack.append(data_subtrack_p)
    if len(Glycosylation) >0:
        data_subtrack.append(data_subtrack_g)
    if len(SUMOylation) >0:
        data_subtrack.append(data_subtrack_su)
    if len(Acetylation) >0:
        data_subtrack.append(data_subtrack_ac)
    if len(Pyrrolidone_carboxylic_acid) >0:
        data_subtrack.append(data_subtrack_pc)
    if len(Palmitoylation) >0:
        data_subtrack.append(data_subtrack_pa)
    if len(Hydroxylation) >0:
        data_subtrack.append(data_subtrack_h)

    track = {
           "label": 'PTM',
           "labelType": 'text',
           "data": data_subtrack
        }
    return track

class MusiteDeep:
    def __init__(self, pdbFile, root_dir):
        self.pdbFile = os.path.abspath( pdbFile)
        self.root_dir = os.path.abspath(root_dir)
        self.FaFi = os.path.join(self.root_dir, 'seqs.fasta')
        self.outDir = os.path.join(self.root_dir, self.__class__.__name__ )
        self.outFi =  os.path.join(self.outDir, 'PTM/musiteDeep_results.txt')
        subprocess.run('mkdir -p ' + self.outDir, stdout = subprocess.PIPE, shell=True)
        self.init_cmd = ' '
        self.track = ''
        
    def run(self):
        self.pdb2fasta(self.pdbFile, self.FaFi)
        cmd = self.init_cmd
        cmd += f"conda run -n MusiteDeep python3 /apps_dk/MusiteDeep_web-master/MusiteDeep/predict_multi_batch.py -input {self.FaFi} -output {self.outDir}/PTM/musiteDeep  -model-prefix '/apps_dk/MusiteDeep_web-master/MusiteDeep/models/Phosphotyrosine;Phosphoserine_Phosphothreonine;Methyllysine;Ubiquitination;SUMOylation'"
        print(cmd)
        subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        self.parse()
        
    def parse(self):
        lines = open(self.outFi).read().strip('\n').split('\n')
        res = []
        for li in lines[2:]:
            pos,residue,ptm,cutoff_ptm = li.split('\t')
            if cutoff_ptm !="None":
                for item in cutoff_ptm.split(';'):
                    ptm_type,score = item.split(':')
                    res.append([int(pos),residue, ptm_type, float(score) ])
        self.track  = musiteDeep2protvista(res)
    
    def pdb2fasta(self, pdbFile, FaFi):
        if not os.path.exists(FaFi):
            cmd2 = f'pdb2fasta {pdbFile} > {FaFi}'
            subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)

def run_MusiteDeep(pdbFile, root_dir):
    mp = MusiteDeep(pdbFile, root_dir)
    mp.run()
    print(mp.track)
    return mp.track
class Main:
    def run(self,pdbFile,root_dir):
        run_MusiteDeep(pdbFile, root_dir)
if __name__ == "__main__":
    import fire
    fire.Fire(Main)