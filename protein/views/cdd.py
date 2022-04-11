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


@cache_page(60 * 15)
def search(request):
    if request.method == "GET":
        this_uuid = request.GET.get('uuid')
        print(this_uuid)
        if this_uuid and cache.get((this_uuid + '_querySet')):
            print("get cache")
            newQ = cache.get((this_uuid + '_querySet'))
            print(newQ)
            protin_id_all = [q['protin_id'] for q in newQ]
            data = {"status": 200}
            data["protin_ids"]=protin_id_all
            return JsonResponse(data)
    else:
        body_unicode = request.body.decode('utf-8')
        body = json.loads(body_unicode)
        print(body)
        newQ = []
        data = {"status": 200}
        reseach = 1
        newParams = json.loads(body_unicode)
        
        # del newParams['pageSize']; del newParams['currentPage']; newParams['field'];  newParams['order']
        this_uuid = body['uuid']
        print(newParams)
        # print(cache.get(this_uuid) )
        if this_uuid and cache.get(this_uuid) and cache.get((this_uuid + '_querySet')):
            print("using cache")
            oldParams = cache.get(this_uuid)
            print(newParams == oldParams )
            # print(oldParams, newParams)
            if newParams == oldParams:
                print("using cache...")
                reseach = 0
                newQ = cache.get((this_uuid + '_querySet'))
        if reseach:
            print('searching ...')
            prot_lenMax = 100000000
            prot_lenMin = int(body['target_len']['min'])
        
            prot_lenMax = int(body['target_len']['max'])
            querySet = protCDncbi.objects.all().filter(length__gt = prot_lenMin, length__lt=prot_lenMax)
            
            domains_keep = slim_domain(body['domains'])
            q_candidate,q_exclude = [],[]
            domain_num =0
            cdWord = []
            for domain in domains_keep:
                keyword =  domain['value'].strip()
                if keyword  and not domain["exclude"]:
                    dCount = domain['count']
                    print("contained key word:--->",keyword)
                    domain_num +=1
                    cd = CDD.objects.filter(cdd_name__icontains=keyword)
                    cdWord.append([domain, cd.count(), cd])
                    # 关键词数, 开发是关掉
                    print("关键词数:", keyword, "CD_number:", cd.count(), cd.values_list('cdd_name', flat=True))
            
            cdWord = sorted(cdWord, key=lambda k: k[1])
           
            # domain_num = 0
            # if len(cdWord) >1:
            for cdw_ in cdWord:
                print("#  round:--->",cdw_)
                domain, keyCount, cd  =  cdw_
                dCount = domain['count']
                q_domain = querySet.filter(cdd_id__in = {q.cdd_id for q in cd})
                # 关键词在protin region 中出现的次数
                # print("first round, contains",  q_domain.count())
                if dCount >1:
                    q_protins = q_domain.values("protin_id").annotate(counts=Count("cdd_name")).filter( counts__gte = dCount )
                    protins = {q["protin_id"] for q in q_protins}
                    print("protins numbers", len(protins))
                else:
                    protins = q_domain.values_list("protin_id",flat=True).distinct()
                    print("protins numbers", len(protins))

                querySet = querySet.filter(protin_id__in = protins)
                print("first round proteins", querySet.count())
           
            
            print("#过滤否定的关键词 filter Out")
            for domain in domains_keep:
                if domain['value']  and domain["exclude"]:
                    dCount = domain['count']
                    keyword =  domain['value'].strip()
                    print("filter out key word:--->",keyword)
                    cd = CDD.objects.filter(cdd_name__icontains=keyword)
                    # 关键词数, 开发是关掉
                    print("# 否定 CD number", cd.count())
                    nr_domian = querySet.filter(cdd_id__in = cd.values_list('cdd_id', flat=True) )
                    querySet = querySet.exclude(protin_id__in =  nr_domian.values_list('protin_id', flat=True).distinct() )
                   
            print("filter ok region:", querySet.count())
            print("filter ok, protin_id:",querySet.values('protin_id').distinct().count())
            # querySet = querySet.exclude(protin_id__in=id_exclude)
            newQ= []
            
            print("query protCDncbiOne")
            newQ = protCDncbiOne.objects.filter(protin_id__in = querySet.values('protin_id').distinct() )
            # save cache
            print("save cache")
            cache.clear()
            this_uuid = str(uuid.uuid4())
            newParams['uuid'] = this_uuid
            myuuid_q = this_uuid + '_querySet'
            cache.set(this_uuid, newParams, 60 * 5)
            cache.set(myuuid_q, newQ, 60 *5)
        
        count= len(newQ)
        print("newQ",count)
        pageSize = 10
        currentPage =1
        order = body['order'];field = body['field']
        requestData = newQ[(currentPage - 1) * pageSize : currentPage * pageSize]
        content = []
        for q in requestData:
            cdd_notes, cdd_names = [],[]
            cdd_ids = q.cdd_idCat.split(', ')
            protin = NrInfo.objects.filter(protin_id = q.protin_id).first()
            for cdd_id_ in cdd_ids:
                qcd = CDD.objects.get(cdd_id = cdd_id_)
                cdd_notes.append(qcd.cdd_desc)
                cdd_names.append(qcd.cdd_name)
            # cdd_notes = '^^'.join(cdd_notes)
            items = {
                "protin_id":q.protin_id,
                "length": protin.length,
                "cdd_names": cdd_names,
                "cdd_notes":cdd_notes,
                "cdd_ids":cdd_ids
            }
            content.append(items)
        
        data["data"] =  content
        data['totalCount'] = count
        data['uuid'] =  this_uuid
        return JsonResponse(data )

def get_all_protin_id(request):
    if request.method == "GET":
        this_uuid = request.GET.get('uuid')
        print(this_uuid)
        if this_uuid and cache.get((this_uuid + '_querySet')):
            print("get cache")
            newQ = cache.get((this_uuid + '_querySet'))
            protin_id_all = [q['protin_id'] for q in newQ]
            data = {"status": 200}
            data["protin_ids"]=protin_id_all
            return JsonResponse(json.dumps(data),safe=False)

def pages(request):
    if request.method == "GET":
        this_uuid = request.GET.get('uuid')
        print(this_uuid)
        if this_uuid and cache.get((this_uuid + '_querySet')):
            print("get cache")
            newQ = cache.get((this_uuid + '_querySet'))
            count = len(newQ)
            print("newQ",len(newQ))
            data = {"status": 200}
            pageSize = int(request.GET.get('pageSize'))
            currentPage = int(request.GET.get('currentPage'))
            order =request.GET.get('order')
            field = request.GET.get('field')
            #  order
            myuuid_order = this_uuid + f'{order}.{field}'
            print("order field", order, field)
            # if myuuid_order !=cache.get(myuuid_order):
            #     if order == "descending":
            #         newQ = sorted(newQ, key=lambda k: k[field], reverse=True)
            #     else:
            #         newQ = sorted(newQ, key=lambda k: k[field])
                
            #     cache.set(myuuid_order, f'{order}.{field}', 60 *60)
            
            requestData = newQ[(currentPage - 1) * pageSize : currentPage * pageSize]
            content = []
            for q in requestData:
                cdd_notes, cdd_names = [],[]
                cdd_ids = q.cdd_idCat.split(', ')
                length = NrInfo.objects.get(id = q.protin_nr_id).length
                for cdd_id_ in cdd_ids:
                    qcd = CDD.objects.get(cdd_id = cdd_id_)
                    cdd_notes.append(qcd.cdd_desc)
                    cdd_names.append(qcd.cdd_name)
                # cdd_notes = '^^'.join(cdd_notes)
                items = {
                    "protin_id":q.protin_id,
                    "length":length,
                    "cdd_names": cdd_names,
                    "cdd_notes":cdd_notes,
                    "cdd_ids":cdd_ids
                }
                content.append(items)
            
            data["data"] =  content
            data['totalCount'] = count
            data['uuid'] =  this_uuid
            return JsonResponse(data )

## 逻辑
def slim_domain(domains):
    new = []
    for dm in domains:
        if dm['value']:
            new.append(dm)
    return new

def filter_cdd_save(param, data):
    target_len_min= float(param['target_len']['min'])
    target_len_max= float(param['target_len']['max'])
    print( target_len_max )
    dt=[]
    print(param)
    field_dt = {
        'protin_id':0,
        "cdd_nameCat":1
    }
    for mydt in data:
        if  mydt['target_len'] >= target_len_min and mydt['target_len'] <= target_len_max:
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
                            # print('match',match.groups())
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
    return dt
