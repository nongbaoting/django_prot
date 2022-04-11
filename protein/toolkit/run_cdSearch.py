from os import path
import re, uuid, pickle
import time
import json
from django.db.models import Count

from collections import defaultdict
from django.utils import timezone
from protein.models import *
from protein.toolkit import myFunctions

cdd_dump_dir = "/dat1/nbt2/proj/21-prot/web/data/res/cdd_dump"

def search_save(params):
    print(params)
    print('searching ...')
    prot_lenMax = 100000000
    prot_lenMin = int(params['target_len']['min'])
    prot_lenMax = int(params['target_len']['max'])

    querySet = protCDncbi.objects.all().filter(length__gt = prot_lenMin, length__lt=prot_lenMax)
    domains_keep = slim_domain(params['domains'])
    domain_num =0
    cdWord = []
    for domain in domains_keep:
        keyword =  domain['value'].strip()
        if keyword  and not domain["exclude"]:
            dCount = domain['count']
            print("contained key word:---> " + keyword)
            domain_num +=1
            cd = CDD.objects.filter(cdd_name__icontains=keyword)
            cdWord.append([domain, cd.count(), cd])
            # 关键词数, 开发是关掉
            print("关键词数: " + keyword + "CD_number: " + str(cd.count()), cd.values_list('cdd_name', flat=True))
    cdWord = sorted(cdWord, key=lambda k: k[1])
    for cdw_ in cdWord:
     
        domain, keyCount, cd  =  cdw_
        dCount = domain['count']
        q_domain = querySet.filter(cdd_id__in = { q.cdd_id for q in cd})
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
            
    # print("filter ok region:", querySet.count())
    # print("filter ok, protin_id:",querySet.values('protin_id').distinct().count())
    # querySet = querySet.exclude(protin_id__in=id_exclude)
    print("query protCDncbiOne")
    newQ = protCDncbiOne.objects.filter(
        protin_id__in = querySet.values('protin_id').distinct() ).values_list(
        'protin_id','cdd_nameCat','cdd_idCat'
    )
    
    tmpFiName = params["uuid"] + '.npz'
    outFi = path.join(cdd_dump_dir, tmpFiName)
    print(outFi)
    myFunctions.numpy_save2file(newQ, outFi)
    return outFi

            
    

## 逻辑
def slim_domain(domains):
    new = []
    for dm in domains:
        if dm['value']:
            new.append(dm)
    return new

