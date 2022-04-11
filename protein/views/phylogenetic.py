
from os import path
import re,uuid
import time
import json
from django.db.models import Count
from protein.models import *
from collections import defaultdict
from django.utils import timezone
from django.core.files import File
from django.http import JsonResponse,FileResponse
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.decorators.cache import cache_page
from django.core.cache import cache
from protein import tasks
from django.http import Http404
from protein.toolkit import run_phylogenetic

# result dir
fasta_dir = "/dat1/nbt2/proj/21-prot/web/data/res/fasta"
phlyo_dir = "/dat1/nbt2/proj/21-prot/web/data/res/phylotree"

# view ---> taks --> utils

def get_fasta(request):
    params = load_POST(request)
    params['work_dir'] = fasta_dir
    seqIDs = params['seqIDs']
    if len(seqIDs) <200:
        myuuid = str(uuid.uuid4())
        params['uuid'] = myuuid
        Fi = run_phylogenetic.getNrFasta(params)
    else:
        Fi = tasks.getNrFasta.apply_async(args = [params])
    if path.exists(Fi):
        fh = open(Fi, 'rb')
        return FileResponse(fh)
    else: 
        raise Http404("request File does not exist! " + Fi)

def run_phylo(request):
    if request.method == "POST":
        params = load_POST(request)
        params['work_dir'] = phlyo_dir
        res = tasks.run_phylo.apply_async(args = [params])
        print(res.task_id, res.status)
        data = {"status": res.status,
            "task_id":res.task_id
        }
        SubmitInfoNew.objects.create(
            uuid = res.task_id,
            proj_type ="Phylogenetic Analysis",
            email_addr = params["email"],
            job_name = params["job_name"],
            upload_date = datetime.datetime.now(),
            tools='',
            params='',
            task_status= res.status
        )
        return JsonResponse(data)
    elif request.method == "GET":
        pass

def tree_files(request):
    if request.method=="GET":
        myuuid =request.GET.get("uuid")
        fileName = ''
        wantedType = request.GET.get("wantedType").strip()
        print(wantedType)
        if wantedType == "tree_pdf":
            fileName = "ggtree.pdf"
        elif wantedType =="msa":
            fileName = "seq.msa"
        elif wantedType == "tree_file":
            fileName = 'seq.msa.iqtree'
        else:
            raise Http404("retrieve wanted type not found!" + wantedType)
        Fi = path.join(phlyo_dir, myuuid, fileName)
        print(Fi)
        if path.exists(Fi):
            print(Fi)
            fh = open(Fi, 'rb')
            return FileResponse(fh)
        else:
            raise Http404("request File does not exist! " + Fi)

## 逻辑
def load_POST(request):
    body_unicode = request.body.decode('utf-8')
    body = json.loads(body_unicode)
    return body

