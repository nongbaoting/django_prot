#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : structure_prediction.py
# @Date            : 2023/07/13 19:04:02
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip,subprocess
import shutil
from django.utils import timezone
from protein.toolkit.modules import PDB

class RoseTTAFold2NA:

    def __init__(self,work_dir,params):
        self.job_name = params['job_name']
        self.program = "RF2NA"
        self.proj_dir = os.path.abspath(work_dir)
        self.work_dir = os.path.join(self.proj_dir, self.program)
        self.outdir   = os.path.join(self.work_dir,  'out')
        os.system('mkdir -p ' + self.outdir)
        self.init_cmd = f'ssh -p3389 nong@172.22.148.150 "source ~/.zshrc;cd {self.work_dir}; conda run -n RF2NA /dat1/apps/RoseTTAFoldNA/RoseTTAFold2NA/run_RF2NA.sh '
    
    def clear(self,):
        cmd = f'ssh -p3389 nong@172.22.148.150 "rm -rf {self.outdir}"'
        subprocess.run(cmd, shell=True)
    
    def format_result(self):
        
        pdb2cif = PDB.Main()
        pdb2cif.scanAF_pdb2cif(self.work_dir)
        # os.system(f"python /dat1/nbt2/pipe/21-prot/bin/PDB.py scanAF_pdb2cif {self.work_dir}")
        ## 如果全变远程,可以在这里cmd加ssh
        
        
        cmd =  f'cd {self.proj_dir}; tar -czvf {self.job_name}.{self.program}.tgz {self.program}/models;'
        cmd += f'cd {self.work_dir}; ln -s models/model_00.cif ranked_0.cif'
        subprocess.run(cmd, shell=True)

    def run(self,protein_fasta, nucleic_type, nucleic_fasta):
        cmd =self.init_cmd + f' {self.outdir} {protein_fasta} {nucleic_type}:{nucleic_fasta}; cp -r {self.outdir}/models ."'
        subprocess.run(cmd, shell=True)
        # shutil.copy(f'{self.outdir}/models', self.work_dir)
        
        self.format_result()
        # self.clear()

class AlphaFold2:
    def __init__(self,work_dir, params):
        self.params = params
        self.job_name = params['job_name']
        self.program = 'AlphaFold2'
        self.proj_dir = os.path.abspath(work_dir)
        self.work_dir = os.path.join(self.proj_dir, self.program)
        self.model_dir =  os.path.join(self.work_dir, 'models')
        self.outdir   = os.path.join(self.work_dir,  'out')
        os.system('mkdir -p ' + self.outdir)
        os.system('mkdir -p ' + self.model_dir)
        self.init_cmd = f'ssh -p3389 nong@172.22.148.150 "source ~/.zshrc;cd {self.work_dir};'
    
    def run(self, faFile):
        status = 1
        params = self.params
        if "AlphaFold2" in params['platform'] or "AlphaFold multimer" in params['platform']:
            model_preset = "monomer"
            if "Multimer" in params['platform']:
                model_preset = 'multimer'
            
            log =self.work_dir + '/run_err.log'
            # cmd = ". /training/nong/app/miniconda3/etc/profile.d/conda.sh; conda activate alphafold; "
            cmd = self.init_cmd + f'conda run -n alphafold python3 /training/nong/protein/apps/alphafold/docker/run_docker.py --gpu_devices 0 --fasta_paths={faFile} --max_template_date=2023-3-9 \
                --data_dir=/training/nong/protein/db/alphafold_dat/  --model_preset={model_preset} --output_dir={self.outdir} 2>{log}"'
            print(f"running alphafold 2,model_preset:{model_preset} =============================>>>>>>>>>>>>>>>>>>>>>>> ", timezone.now())
            run = subprocess.run(cmd, shell=True)
           
            pdb2cif = PDB.Main()
            pdb2cif.scanAF_pdb2cif(self.outdir)
            print("alphafold returncode: >>>>>>>>>>>>>>> ", run.returncode)

            if run.returncode == 0:
                if os.path.exists(f'{self.outdir}/protein/unrelaxed_model_1_multimer.cif'):
                    os.system(f'cp {self.outdir}/protein/unrelaxed_model_1_multimer.cif  {self.model_dir}/unrelaxed_model_1.cif')
                os.system(f"cp -r {self.outdir}/protein/rank*   {self.model_dir}/")
                os.system(f"cp -r {self.outdir}/protein/unrelaxed_model_1.cif   {self.model_dir}/unrelaxed_model_1.cif")
                cmd =  f'cd {self.proj_dir}; tar -czvf {self.job_name}.{self.program}.tgz {self.program}/models;'
                cmd += f'cd {self.work_dir}; ln -s models/ranked_0.cif ranked_0.cif; ln -s models/unrelaxed_model_1.cif unrelaxed_model_1.cif'
                os.system(cmd)
                
            status = run.returncode
        return status

class Main:
    def RF2NA(self, work_dir,protein_fasta, nucleic_type, nucleic_fasta):
        rf2na = RoseTTAFold2NA(work_dir)
        rf2na.run(protein_fasta, nucleic_type, nucleic_fasta)

if __name__ == '__main__':
    fire.Fire(Main)