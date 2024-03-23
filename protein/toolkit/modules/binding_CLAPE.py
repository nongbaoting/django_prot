import os, subprocess, fire


def parse2protvista(outFi,label):
    fragments = []
    lines = open(outFi, 'r').read().split('\n')
    seqs = list(lines[1])
    sites = list(lines[2])
   
    for index, char in enumerate(sites):
        if char=='1':
            fragment_ = {
                "start": index + 1,
                "end": index + 1, 
                "tooltipContent": f"Position: {index}, Seq: {seqs[index]}", 
                "color": '#DC143C' 
            }
            fragments.append(fragment_)
    
    fragments_dt = {"fragments":fragments}
    data_subtrack = {
                    "accession": 'CLAPE' ,
                    "type": f'site',
                    "label": "CLAPE" ,
                    "shape": 'circle',
                    "labelTooltip": f"CLAPE {label} binding sites" ,
                    "locations": [fragments_dt]
                    }
    
    track = {
           "label": f'{label} binding sites',
           "labelType": 'text',
           "data":[data_subtrack]
        }
    
    return track

class CLAPE:
    def __init__(self,pdbFile, root_dir):
        self.init_cmd = ''
        self.root_dir = os.path.abspath(root_dir)
        self.pdbFile = os.path.abspath( pdbFile)
        self.name = os.path.basename(pdbFile).split('.')[0]
        self.FaFi = os.path.join(self.root_dir,'seqs.fasta')
        self.out_dna = os.path.join(self.root_dir, 'CLAPE_DNA.txt')
        self.out_rna = os.path.join(self.root_dir, 'CLAPE_RNA.txt')
        self.out_ab  = os.path.join(self.root_dir, 'CLAPE_AB.txt')

    def run_clape(self):
        self.pdb2fasta()
        cmd = self.init_cmd
        cmd += f"conda run -n CLAPE python /apps_dk/CLAPE/clape.py  --ligand DNA --input {self.FaFi} --output {self.out_dna};"
        cmd += f"conda run -n CLAPE python /apps_dk/CLAPE/clape.py  --ligand RNA --input {self.FaFi} --output {self.out_rna};"
        cmd += f"conda run -n CLAPE python /apps_dk/CLAPE/clape.py  --ligand AB  --input {self.FaFi} --output {self.out_ab}"
       
        print(cmd)
        subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        return self.parse_clape() 
    def parse_clape(self):
        res_dna = parse2protvista(self.out_dna, "DNA")
        res_rna = parse2protvista(self.out_rna, "RNA")
        res_ab = parse2protvista(self.out_ab, "Antibody")
        return [res_dna, res_rna,res_ab ]
    def pdb2fasta(self):
        if not os.path.exists(self.FaFi):
            cmd2 = f'pdb2fasta {self.pdbFile} > {self.FaFi}'
            subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)
def run_clape(pdbFile,root_dir):
    clape = CLAPE(pdbFile, root_dir)
    track = clape.run_clape()
    return track

class Main:
    def run_clape(self,pdbFile,root_dir):
        import json
        clape = CLAPE(pdbFile, root_dir)
        track = clape.run_clape()
        with open('track.json', 'w') as fp:
            json.dump(json.dumps(track[0]), fp)

if __name__ == '__main__':
    print('')
    fire.Fire(Main)