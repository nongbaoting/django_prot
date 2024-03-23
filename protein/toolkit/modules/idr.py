
import fire, os, subprocess, re, json
import pandas as pd

def pdb2fasta(pdbFile, FaFi):
    if not os.path.exists(FaFi):
        cmd2 = f'pdb2fasta {pdbFile} > {FaFi}'
        subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)

def flDPnn2protvista(outFi):
    sites,bins,scores = [],[],[]
    with open(outFi,'r') as f:
        for li in f:
            if re.match('>', li):
                id = li.strip().split('\s')[0][1:]
                sites = next(f).strip().split(',')
                bins = next(f).strip().split(',')
                scores = next(f).strip().split(',')
                break
    fragments = []
    print(sites)

    start, end = 0,0
    fragment_score = []
    for index, char in enumerate(bins):
        # print(char)
        if char == "1" and end ==0:
            fragment_score.append(scores[index])
            start = index +1
            end = 1
        elif char =="0" and start !=0:
            end = index 
            fragment_ = {
                    "start": start,
                    "end": end, 
                    "tooltipContent": f"IDR <br>Range: {start} - {end}", 
                    "color": '#BC80BD', 
                    "score": fragment_score
                    }
            start,end = 0,0
            fragment_score = []
            fragments.append(fragment_)
    
    fragments_dt = {"fragments":fragments}
    data_subtrack = {
                    "accession": 'IDR' ,
                    "type": 'Domain',
                    "label": "flDPnn",
                    "labelTooltip":"IDR",
                    "locations": [
                        {"fragments":fragments}             
                        ]
                    }

    track = {
        "label": 'IDR',
        "labelType": 'text',
        "data":[data_subtrack]
    }

    return track


class flDPnn:
    def __init__(self, pdbFile, root_dir):
        self.pdbFile = os.path.abspath( pdbFile)
        self.root_dir = os.path.abspath(root_dir)
        self.FaFi = os.path.join(self.root_dir, 'seqs.fasta')
        self.outDir = os.path.join(self.root_dir, self.__class__.__name__ )
        self.outFi =  os.path.join(self.root_dir, self.__class__.__name__, 'results.csv')
        subprocess.run('mkdir -p ' + self.outDir, stdout=subprocess.PIPE, shell=True)
        self.init_cmd = f'cd {self.outDir}; ' 
        self.res = []
        
    def run(self):
        pdb2fasta(self.pdbFile, self.FaFi)
        cmd = f"cd {self.outDir}; conda run -n flDPnn python /apps_dk/fldpnn-master/run_flDPnn.py {self.FaFi}"
        # subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        self.parse()
        
    def parse(self):
        self.res =  flDPnn2protvista(self.outFi)
    
    def save(self,outJson):
        with open(outJson, 'w') as fp:
            json.dump(json.dumps(self.res), fp)

def run_idr(pdbFile, root_dir):
    lms = flDPnn(pdbFile, root_dir)
    lms.run()
    print(lms.res)
    return lms.res

class Main:
    def run(self,fasta):
        from mymetal import mbp,iof
        metal = Mymetal(fasta)
        metal.run()
        metal.save
    def run_idr(self,pdbFile, root_dir):
        lms = flDPnn(pdbFile, root_dir)
        lms.run()
        lms.save("flDPnn.json")

if __name__ == "__main__":
    fire.Fire(Main)