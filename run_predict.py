#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : run_predict.py
# @Date            : 2021/08/08 13:03:37
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import django
import os
import sys
import fire
import subprocess
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import datetime
import json
from collections import defaultdict
from django.utils import timezone

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'django_prot.settings')
if django.VERSION >= (1, 7):
    django.setup()


class MyHandler(FileSystemEventHandler):
    # def on_modified(self, event):
    #     print("文件被修改了 %s" % event.src_path)
    def on_created(self, event):
        print("文件被创建了 %s" % event.src_path)
        run = RUN()
        if os.path.isfile(event.src_path):
            # running_date = datetime.datetime.now()
            run.run_struct_predict(faFile=event.src_path)
            # completed_date = datetime.datetime.now()
            # run.stat(run.pre, running_date, completed_date)


def myconverter(o):
    if isinstance(o, datetime.datetime):
        return o.__str__()


uploads_47 = "/training/nong/web/data/uploads/"
web_base = '/training/nong/protein/db/web/outputs/results/'
res_base = '/training/nong/protein/res/'


class RUN:
    def stat(self, job_name, run_date, completed_date, statFi="/training/nong/web/prot/protein/static/stats/structure.json"):
        dat = defaultdict(dict)
        if os.path.getsize(statFi):
            with open(statFi, 'r') as f:
                dat = json.load(f)

        if job_name not in dat:
            dat[job_name] = {'job_name': job_name,
                             "running_date":  run_date,
                             'completed_date':  completed_date,
                             }

        with open(statFi, 'w') as fo:
            json.dump(dat, fo, default=myconverter)
        print(dat)

    def alphafold2(self, faFile, params):
        if "AlphaFold 2" in params['platform'] or "AlphaFold multimer" in params['platform']:
            model_preset = "monomer"
            if "Multimer" in params['platform']:
                model_preset = 'multimer'
            pre = os.path.basename(faFile).split('.fa')[0]
            res_root = res_base + 'alphafold/'
            res_dir = res_root + pre
            os.system('mkdir -p ' + res_dir)
            web_dir = web_base + 'alphafold/' + pre
            log = res_dir + '/run_err.log'
            cmd = ". /training/nong/app/miniconda3/etc/profile.d/conda.sh; conda activate alphafold; "
            cmd += f"python3 /training/nong/protein/apps/alphafold/docker/run_docker.py --fasta_paths={faFile} --max_template_date=2022-3-9 \
                --data_dir=/training/nong/protein/db/alphafold_dat/  --model_preset={model_preset} --output_dir={res_root} 2>{log}; "
            print(f"running alphafold 2,model_preset:{model_preset} =============================>>>>>>>>>>>>>>>>>>>>>>> ", timezone.now())
            run = subprocess.run(cmd, shell=True)
            print("alphafold returncode: >>>>>>>>>>>>>>> ", run.returncode)

            if run.returncode == 0:
                if os.path.exists(f'{res_dir}/unrelaxed_model_1_multimer.pdb'):
                    os.system(f'cp {res_dir}/unrelaxed_model_1_multimer.pdb  {res_dir}/unrelaxed_model_1.pdb')
                os.system(
                    f'cd {res_root}; tar -czvf {pre}_alphafold.tar.gz {pre}/rank*pdb {pre}/unrelaxed_model_1.pdb')
                os.system(f'mkdir -p {web_dir}/{pre}')
                os.system(f'cp {res_root}/{pre}_alphafold.tar.gz {web_dir}')

                os.system(
                    f"cp -r {res_dir}/rank*pdb  {res_dir}/unrelaxed_model_1.pdb {web_dir}/{pre}")
            self.alphafold2_code = run.returncode
        else:
            self.alphafold2_code = 1

    def roseTTAFold(self, faFile, params):
        if "RoseTTAFold" in params['platform']:
            roseMode = params['RoseTTAFold_mode'] 
            bins = {
                'pyrosetta': 'run_pyrosetta_ver.sh',
                'e2e': 'run_e2e_ver.sh',
            }
            cmd = ". /training/nong/app/miniconda3/etc/profile.d/conda.sh; conda activate RoseTTAFold;"
            pre = os.path.basename(faFile).split('.fa')[0]
            res_root = res_base + "roseTTAFold/"
            res_dir = res_root + pre
            web_dir = web_base + 'roseTTAFold/' + pre
            log = res_dir + '/run_err.log'
            os.system(f'mkdir -p ' + res_dir)
            cmd += "/training/nong/protein/apps/RoseTTAFold/" + \
                bins[mode] + ' ' + faFile + ' ' + res_dir + ' 2>' + log

            print("running RoseTTAFold= =============================>>>>>>>>>>>>>>>>>>>>>>> ", timezone.now())
            run = subprocess.run(cmd, shell=True)

            cmd_res = f'cd {res_dir}; mkdir -p model_rel; cp model/* model_rel/;rm -rf model; mv model_rel model;'

            print("RoseTTAFold returncode: >>>>>>>>>>>>>>> ", run.returncode)
            if run.returncode == 0:
                os.system(f'mkdir -p {web_dir}/{pre}')
                if bins[mode] == 'run_pyrosetta_ver.sh':
                    os.system(cmd_res)
                    os.system(
                        f'cd {res_root}; tar -czvf {pre}_roseTTAFold.tar.gz {pre}/model')
                    os.system(f'cp -r {res_dir}/model {web_dir}/{pre}/')
                else:
                    # end to end
                    os.system(
                        f'cd {res_root}; tar -czvf {pre}_roseTTAFold.tar.gz {pre}/t000_.e2e.pdb')
                    os.system(f'cp -r {res_dir}/t000_.e2e.pdb {web_dir}/{pre}/')

                os.system(f'cp -r  {res_root}/{pre}_roseTTAFold.tar.gz {web_dir}')

            self.roseTTAFold_code = run.returncode
        else:
            self.roseTTAFold_code =1

    def run_struct_predict(self, faFile):
        from protein.models import SubmitInfo
        self.fa = faFile
        self.basename = os.path.basename(faFile)
        self.pre = self.basename.split('.fa')[0]
        obj = SubmitInfo.objects.get(job_name=self.pre)

        status = ['1', '1']
        print(obj.running_date)

        if obj.running_date is not None and '0' in status:
            print('project have run before, skipping....')
        else:
            obj.running_date = timezone.now()
            obj.save()
            roseMode = 'pyrosetta'
            if obj.params is not None:
                print('running.....')
                params = json.loads(obj.params)
                self.alphafold2(faFile,params)   
                self.roseTTAFold(faFile, params)
                status = [str(self.alphafold2_code), str(self.roseTTAFold_code)]
                
            django.db.close_old_connections()
            obj = SubmitInfo.objects.filter(job_name=self.pre).first()
            obj.run_status = ','.join(status)
            obj.completed_date = timezone.now()
            obj.save()
        django.db.close_old_connections()

    def watchdog(self, mypath=uploads_47 + 'structure_predict'):
        event_handler = MyHandler()
        observer = Observer()
        observer.schedule(event_handler, mypath, recursive=True)
        observer.start()
        try:
            while True:
                # every ten seconds
                time.sleep(10)
                # print('hello')
        except:
            observer.stop()
        observer.join()

    def watch_blast(self):
        from myscripts import run_blast
        run_blast.watch_dog()


if __name__ == "__main__":

    fire.Fire(RUN)
