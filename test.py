
# Create your tests here.
import django
import os, re
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'django_prot.settings')
django.setup()
from protein.toolkit import *
import json, fire
from multiprocessing import Pool

def scanAndFind_pattern(mydir, mypattern):
    wantFiles = []
    for entry in os.scandir(mydir):
        if (entry.is_file() or entry.is_link()) and mypattern.search(entry.name):
            wantFiles.append(entry)
        elif entry.is_dir():
            wantFiles.extend( scanAndFind_pattern(entry.path, mypattern) )
    return wantFiles

class Main:

    def test1(self,paraFile):
        '''
        python ./test.py test1 /apps/resultData/pdb_annotations/results/430eaaa9-400e-428c-92b9-c388b70a5453/params.json >test.log 2>&1 &
        '''
        with open(paraFile) as f:
            params = json.loads(json.load(f))
            print(params)
            results = run_pdb_domain_annotations.domain_annotations( params )
            
    def run_pdb_annotation(self, pdbDir, outDir,core=12):
        re_pdb = re.compile('.pdb$')
        pool = Pool(core)
        res = []
        for entry in scanAndFind_pattern(pdbDir,re_pdb):
            name = entry.name.split('-')[1]
            params ={}
            params['job_name'] = name
            params['chain'] = 'A'
            params['uuid'] = name
            params['proj_type'] = "PDB Domain Annotations"
            params['work_dir'] = os.path.join(outDir, 'results')
            params['pdb_annotations_dir_remote'] = ''
            params['pdbfile'] =  entry.path
            params['ip'] = '172.22.148.191'
            res_ = pool.apply_async(run_pdb_domain_annotations.domain_annotations, (params,))
        pool.close()
        pool.join()
        
if __name__ == '__main__':
    fire.Fire(Main)