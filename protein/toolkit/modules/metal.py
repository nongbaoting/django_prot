
import fire, os,subprocess
import pandas as pd

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
    print(sites)
    for index, char in enumerate(sites):
        if char == "1":
            fragment_ = {
                "start": index + 1,
                "end": index + 1, 
                "tooltipContent": f"Position: {index}, Probability: {scores[index]}", 
                "color": '#DC143C'
            }
            fragments.append(fragment_)
    
    fragments_dt = {"fragments":fragments}
    data_subtrack = {
                    "accession": 'LMetalSite' ,
                    "type": 'site',
                    "label": "LMetalSite" ,
                    "shape": 'circle',
                   
                    "locations": [fragments_dt]
                    }
    
    track = {
           "label": f'{label}',
           "labelType": 'text',
           "data":[data_subtrack]
        }
    
    return track


class LMetalSite:
    def __init__(self, pdbFile, root_dir):
        self.pdbFile = os.path.abspath( pdbFile)
        self.root_dir = os.path.abspath(root_dir)
        self.FaFi = os.path.join(self.root_dir, 'seqs.fasta')
        self.outDir = os.path.join(self.root_dir, self.__class__.__name__ )
        self.outFi =  os.path.join(self.root_dir, self.__class__.__name__, 'seqs_predictions.csv')
        subprocess.run('mkdir -p ' + self.outDir, stdout=subprocess.PIPE, shell=True)
        self.init_cmd = 'ssh -p3389 nong@172.22.148.150 "source ~/.zshrc;'
        self.res = []
    def run(self):
        pdb2fasta(self.pdbFile, self.FaFi)
        cmd = self.init_cmd
        cmd += f"conda run -n LMetalSite python /dat1/apps/LMetalSite/script/LMetalSite_predict.py --fasta {self.FaFi} --outpath {self.outDir}"
        cmd +='"'
        print(cmd)
        subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        self.parse()
    def parse(self):
        dt = pd.read_csv(self.outFi,dtype=str)
        print(dt.columns)
        ZN = lMetalSite2protvista(dt.ZN_prob[0], dt.ZN_pred[0],'Zn')
        CA = lMetalSite2protvista(dt.CA_prob[0], dt.CA_pred[0],'Ca')
        MG = lMetalSite2protvista(dt.MG_prob[0], dt.MG_pred[0],'Mg+')
        MN = lMetalSite2protvista(dt.MN_prob[0], dt.MN_pred[0],'Mn+')
        self.res =  [ZN,CA,MG,MN]
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