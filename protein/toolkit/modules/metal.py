
import fire,os,subprocess

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

class LMetalSite:
    def __init__(self, pdbFile, root_dir):
        self.pdbFile = os.path.abspath( pdbFile)
        self.root_dir = os.path.abspath(root_dir)
        self.FaFi = os.path.join(self.root_dir, 'seqs.fasta')
        self.outDir = os.path.join(self.root_dir, self.__class__.__name__)
        self.outFi =  os.path.join(self.root_dir, self.__class__.__name__, 'seqs_predictions.csv')
        subprocess.run('mkdir -p ' + self.outDir, stdout=subprocess.PIPE, shell=True)
        self.init_cmd = 'ssh -p3389 nong@172.22.148.150 "source ~/.zshrc;'
    
    def run(self):
        pdb2fasta(self.pdbFile, self.FaFi)
        cmd = self.init_cmd
        cmd += f"conda run -n LMetalSite python /dat1/apps/LMetalSite/script/LMetalSite_predict.py --fasta {self.FaFi} --outpath {self.outDir}"
        cmd +='"'
        print(cmd)
        subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
    
    def parse(self):
        tab = []
        
def run_LMetalSite(pdbFile, root_dir):
    lms = LMetalSite(pdbFile, root_dir)
    lms.run()

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