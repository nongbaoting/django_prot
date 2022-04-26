from os import path
import re,uuid, pickle
import time
import json
from django.db.models import Count
from collections import defaultdict
from django.utils import timezone
from django.core.files import File
from django.http import JsonResponse, HttpResponse,Http404
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.decorators.cache import cache_page
from django.core.cache import cache
from protein.models import *
from protein.toolkit import run_cdSearch,myFunctions
from protein import tasks
import numpy as np

cdd_dump_dir = "/dat1/nbt2/proj/21-prot/web/data/res/cdd_dump"

def search_save(request):
    # tasks.run_phylo.apply_async(args = [params])
    params = myFunctions.load_POST(request)
    params['proj_type'] = "CD Search"
    params["email"] =''
   
    params["tools"] =''
    res = tasks.cd_search_save.apply_async(args = [params])
    params['uuid']  = res.task_id
    myFunctions.create_submit_form(params, res)
    return JsonResponse({'msg':200})

saveCddDict = defaultdict(dict)
saveCddLt  = []
def retrieve_save(request):
    if request.method == "GET":
        this_uuid = request.GET.get('uuid')
        dumpFi = path.join(cdd_dump_dir, this_uuid + '.npz')
        newQ = load_cdd_data(this_uuid)
        data = {}
        data["data"] =  []
        data['totalCount'] = 0
        if len(newQ) > 0:
            count = len(newQ)
            data = {"status": 200}
            pageSize = int(request.GET.get('pageSize'))
            currentPage = int(request.GET.get('currentPage'))
            order = request.GET.get('order')
            field = request.GET.get('field')
            isSort = request.GET.get('isSort')
            anno_source = request.GET.get('anno_source')
            if isSort=="Yes":
                print("order field", order, field)
                field_dt = {
                    'protin_id':0,
                    "length":2
                }
                field_num = field_dt[field]
                if order == "descending":
                    newQ = sorted(newQ, key=lambda k: k[field_num], reverse=True)
                else:
                    newQ = sorted(newQ, key=lambda k: k[field_num])
                saveCddDict[this_uuid] = newQ

            content = get_one_page(newQ,currentPage, pageSize,anno_source)
            data["data"] =  content
            data['totalCount'] = count
            data['uuid'] =  this_uuid
            return JsonResponse(data )
        else:
            data['uuid'] =  this_uuid
            return JsonResponse(data )

def filter_cdd_save(request):
    param = myFunctions.load_POST(request)
    print(param)
    newQ = load_cdd_data(param['uuid'])
    target_len_min= float(param['target_len']['min'])
    target_len_max= float(param['target_len']['max'])
    anno_source = param['anno_source']
    print( target_len_max )
    
    print(param)
    field_dt = {
        'protin_id':0,
        "length":2
    }
    dt=[]
    for mydt in newQ:
        # if  mydt['target_len'] >= target_len_min and mydt['target_len'] <= target_len_max:
        if True:
            checkDomain =[True]
            for domain in param['domains']:
                if domain['value']:
                    isInclude = False
                    field_num = field_dt[ domain['type'] ]
                    dCount = domain['count']
                    keyword =  domain['value'].strip()
                    if not domain["exclude"]:
                        # print("key word", keyword, "count: ", dCount,mydt[ domain['type']] )
                        ire_key = f"({keyword}.*?)" +  "{" + str(dCount) + ",}"
                        re_k = re.compile(ire_key, re.I)
                        match = re_k.search(mydt[field_num])
                        if match:
                            isInclude = True
                    else:
                        # print('aa')
                        re_exclude = re.compile(keyword, re.I)
                        if not re_exclude.search(mydt[ field_num ]):
                            isInclude = True
                    checkDomain.append(isInclude)
            # print(checkDomain)
            if not False in checkDomain:
                dt.append(mydt) 
    this_uuid = str(uuid.uuid4())
    saveCddDict[this_uuid] = dt
    content = get_one_page(dt,1, 10,anno_source)
    data = {}
    data["data"] =  content
    data['totalCount'] = len(dt)
    data['uuid'] =  this_uuid
    return JsonResponse(data )

## 逻辑
def slim_domain(domains):
    new = []
    for dm in domains:
        if dm['value']:
            new.append(dm)
    return new

def load_cdd_data(this_uuid):
    dumpFi = path.join(cdd_dump_dir, this_uuid + '.npz')
    newQ = []
    if this_uuid in saveCddDict:
        print('using cache...')
        newQ = saveCddDict[this_uuid]
    elif  path.exists(dumpFi):
        print(dumpFi, 'ok')
        newQ = np.load(dumpFi)['cdd']
        print("newQ",len(newQ))
        saveCddLt.append(this_uuid)
        saveCddDict[this_uuid] = newQ
        ## 保持dict keys 数量
        # blastDict, blastLt = shift_dict(blastDict, res_id, blastLt, keep_len = 10)
        if len(saveCddLt) > 10:
            headKey = saveCddLt.pop(0)
            print("delte key:", headKey )
            del saveCddDict[headKey]
            print("keys number: ", len( saveCddDict.keys() ) )
    return newQ

def get_one_page(newQ, currentPage,pageSize,anno_source):
    requestData = newQ[(currentPage - 1) * pageSize : currentPage * pageSize]
    content = []
    # cddall = CDD.objects.all()
    for q in requestData:
        print(q)
        protin_id = q[0]
        cdd_locs, cdd_ids,cdd_notes, cdd_names = [],[],[],[]
        protin = NrInfo.objects.filter(protin_id = protin_id).first()
        regions = []
        if anno_source == "all":
            regions = protCD.objects.filter(protin_id = protin_id)
        elif anno_source == "ncbi":
            regions = protCDncbi.objects.filter(protin_id = protin_id)
        for qr in regions:
            qcd = CDD.objects.get(cdd_id = qr.cdd_id_id)
            cdd_ids.append( qr.cdd_id_id)
            cdd_notes.append(qcd.cdd_desc)
            cdd_names.append(qcd.cdd_name)
            cdd_locs.append([qr.start, qr.end, qcd.cdd_name, qcd.cdd_desc])
        # cdd_notes = '^^'.join(cdd_notes)
        items = {
            "protin_id" :protin_id,
            "length": protin.length, # 改成q[2]
            'desc':protin.desc, # 改成q[1]
            "cdd_names": cdd_names,
            "cdd_notes":cdd_notes,
            "cdd_ids":cdd_ids,
            "cdd_locs":cdd_locs
        }
        content.append(items)

    print("addCDcorlor")
    content = myFunctions.addCDcorlor(content)
    return  content
