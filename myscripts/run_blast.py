#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : run_blast.py
# @Date            : 2022/03/05 14:54:10
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, json, re
from os import path
import time
import datetime
import json

from collections import defaultdict
import django
from django.utils import timezone
from myscripts import parseBLAST_HMM
import subprocess
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

dbDir = {
"blast":{
      'nr': '/dat1/dat/db/nr/all/nr',
      'nr_gt400': '/dat1/dat/db/nr/nr400/nr_gt400',
      'meta': "/dat1/dat/db/nr/meta/db/meta",
 },
"jackhmmer":{
    'nr': "/dat1/dat/db/nr/nr.fa",
    'nr_gt400':'/dat1/dat/db/nr/nr400/nr400.fa',
    'meta': "/dat1/dat/db/nr/meta/meta_genome_prot.fa"
}
}

# os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'django_prot.settings')
# django.setup()

def check_file(Fi):
    filesize = os.path.getsize(Fi)
    if filesize == 0:
        time.sleep(1)
    while True:
        if(filesize == os.path.getsize(Fi)):
            break
        time.sleep(1)
        filesize = os.path.getsize(Fi)
    
        
class MyHandler(FileSystemEventHandler):
    def params(self, regex):
        self.regex = regex
    def on_created(self, event):
        print("file was created: %s" % event.src_path)
        
        if os.path.isfile(event.src_path) and self.regex.search(os.path.basename(event.src_path)):
            check_file( event.src_path)
            print('run: ' + event.src_path)
            # 创建新实例
            run = BLAST(event.src_path)
            run.run()
        
def myconverter(o):
    if isinstance(o, datetime.datetime):
        return o.__str__()

# watch dog include to run predict 
blast_indir = "/dat1/nbt2/proj/21-prot/web/data/uploads/blast"
blast_outdir ="/dat1/nbt2/proj/21-prot/web/data/res/blast/"
def watch_dog(mypath=blast_indir):
    regex = re.compile('.json$')
    event_handler = MyHandler()
    event_handler.params( regex)
    observer = Observer()
    observer.schedule(event_handler, mypath, recursive=True)
    observer.start()
    try:
        while True:
            time.sleep(10)
    except:
        observer.stop()
    observer.join()


class BLAST:
    def __init__(self, in_json):
        self.file_id = os.path.basename(in_json).split('.')[0]
        self.fa_file = os.path.join(os.path.dirname(in_json), f"{self.file_id}.fa")
        self.out_dir = path.join(blast_outdir, self.file_id)
        os.system("mkdir -p "+ self.out_dir)
        f = open(in_json)
        dt = json.load(f)
        # print(in_json)
        self.jackhmmer_params = f"-N { dt['jackhmmer_iter']} " + dt['jackhmmer_advance']
        self.psiblast_params = f"-num_iterations { dt['psiblast_iter']} " + dt['psiblast_advance']
        self.blastp_params = ''
        self.dbName = dt['db']
        self.programs = dt['program']
        f.close()
 
    def run(self):
        from protein.models import SubmitInfoNew
        obj =  SubmitInfoNew.objects.filter(uuid=self.file_id).first()
        if obj and obj.completed_date is not None:
            print('project have run before, skipping....')
        elif obj:
            obj.running_date = timezone.now()
            obj.save()
            status_psi = self.psiblast()
            status_jack =  self.jackhmmer()  
            status_blastp = self.blastp()
            
            django.db.close_old_connections()
            obj =  SubmitInfoNew.objects.filter(uuid=self.file_id).first()
            obj.completed_date = timezone.now()
            obj.run_status = ','.join([str(i)  for i in [status_psi ,status_jack,status_blastp]])
            obj.save()
        else:
            print("DB not found uuid:", self.file_id)
            print( obj)
        django.db.close_old_connections()
        
    def blastp(self):
        program="BLASTP"
        status = 1
        if program in self.programs:
            self.blastp_json = path.join(self.out_dir, 'BLASTP.json' )
            self.blastp_archi_json = path.join(self.out_dir, 'BLASTP_archi.json' )
            self.blastp_out = path.join(self.out_dir,  'BLASTP.xml' )
            cmd=f"blastp -db  {dbDir['blast'][self.dbName]} -query  {self.fa_file} -max_target_seqs 100000 \
            -outfmt 5  -out { self.blastp_out} -num_threads 16 \
             {self.blastp_params} "
            
            print(f"start {program}: >>>>>>>>>>>>>>> \n\n", cmd)
            run = subprocess.run(cmd, shell=True)
            print(f"{program} returncode: >>>>>>>>>>>>>>> ", run.returncode)
            status = run.returncode
            if run.returncode == 0:
               print(f'{program}  success!')
               print(f'{program}  success!')
               ## /parse file
               blast = parseBLAST_HMM.BLAST()
               blast.parse_psiblast(self.blastp_out, self.blastp_json,self.blastp_archi_json )              
        else:
            status = -1
            print("Did not choose " + program)
        return status
        
    def psiblast(self):
        program = "PSI-BLAST"
        status = 1
        if program in self.programs:
            self.psi_out = path.join(self.out_dir,  'psiblast.xml' )
            self.psi_json = path.join(self.out_dir, 'psiblast.json' )
            self.psi_out_archi = path.join(self.out_dir,  'psiblast_archi.json' )
            cmd=f"psiblast -db  {dbDir['blast'][self.dbName]} -query  {self.fa_file} -max_target_seqs 1000000 \
            -outfmt 5  -out { self.psi_out} -num_threads 16 \
             {self.psiblast_params} "
            
            print(f"start {program}: >>>>>>>>>>>>>>> \n\n", cmd)
            run = subprocess.run(cmd, shell=True)
            print(f"{program} returncode: >>>>>>>>>>>>>>> ", run.returncode)
            status = run.returncode
            if run.returncode == 0:
               print(f'{program}  success!')
               print(f'{program}  success!')
               ## /parse file
               blast = parseBLAST_HMM.BLAST()
               blast.parse_psiblast(self.psi_out, self.psi_json ,self.psi_out_archi)              
        else:
            print("Did not choose " + program)
            status = -1
        return status
    
    def jackhmmer(self, ):
        program = "jackhmmer"
        status = 1
        if program in self.programs:
            out_hmm = path.join(self.out_dir ,'out.hmm')
            out_hmmtbl = path.join(self.out_dir ,'out.hmm.tbl')
            out_hmmdomttbl = path.join(self.out_dir ,'out.hmm.domttbl')
            out_json =  path.join(self.out_dir,"jackhmmer.json")
            out_architectures =  path.join(self.out_dir,"jackhmmer_archi.json")
            cmd=f"jackhmmer  {self.jackhmmer_params} \
                -o {out_hmm} --tblout {out_hmmtbl} --domtblout {out_hmmdomttbl} \
                --cpu 16 {self.fa_file} {dbDir['jackhmmer'][self.dbName]}"
            print(f"start {program}: >>>>>>>>>>>>>>> \n\n", cmd)
            run = subprocess.run(cmd, shell=True)
            print(f"{program} returncode: >>>>>>>>>>>>>>> ", run.returncode)
            status = run.returncode
            if run.returncode == 0:
               print(f'{program}  success!')
               ## /parse file
               hmmer = parseBLAST_HMM.HMMER()
               hmmer.parse_jackhmmer(out_hmm, out_hmmdomttbl, out_json,out_architectures ) 
        else:
            print("Did not choose " + program)
            status = -1
        return status
    
if __name__ == '__main__':
    fire.Fire(BLAST)
    