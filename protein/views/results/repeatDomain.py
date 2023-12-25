from django.http import JsonResponse, HttpResponse, FileResponse, Http404

from protein.models import RepeatDomain
from protein.models import *
import os, json
from protein.toolkit import run_cdSearch, myFunctions
from django.views.decorators.cache import cache_page

# @cache_page(60 * 15)
def fetchData(request):

    pageSize = int(request.GET.get("pageSize"))
    currentPage = int(request.GET.get("currentPage"))
    order = request.GET.get("order")
    field = request.GET.get("field")
    
    params =  request.GET.dict()
    print(params)
    print(params['domains'])
    domains =  json.loads( params['domains'])['data']
    info =  RepeatDomain.objects.all()
    for domain in domains:
        print(domain)
        keyword = domain['value'].strip()
        if keyword and not domain["exclude"]:
            dCount = domain['count']
            print("contained key word:--->", keyword)
            info = info.filter(dm_names__icontains=keyword)

    info = info.order_by('min_len')
    totalCount = info.count() 
    requestData = info[(currentPage-1) * pageSize : currentPage * pageSize]
    content = []
    for q in requestData:
        protin_id = q.reprProtein
        cdd_locs, cdd_ids, cdd_notes, cdd_names = [], [], [], []
        protin = NrInfo.objects.filter(protin_id=protin_id).first()
        regions = protCD.objects.filter(protin_id=protin_id)
        for qr in regions:
            qcd = CDD.objects.get(cdd_id=qr.cdd_id_id)
            cdd_ids.append(qr.cdd_id_id)
            cdd_notes.append(qcd.cdd_desc)
            cdd_names.append(qcd.cdd_name)
            cdd_locs.append(
                [qr.start, qr.end, qcd.cdd_name, qcd.cdd_desc])
        # cdd_notes = '^^'.join(cdd_notes)
        items = {
            "q_id":q.id,
            "protin_id": protin_id,
            "length": protin.length,
            "cdd_names": cdd_names,
            "cdd_notes": cdd_notes,
            "cdd_ids": cdd_ids,
            "cdd_locs": cdd_locs,
            "clusters": q.clusters,
            "clusters_num": len(q.clusters.split(":")),
            "numProtein"  : q.numProtein,
            "min_len" : q.min_len,
            "max_len" : q.max_len,
        }
        content.append(items)
    
    data = {"status": 200}
    
    content = myFunctions.addCDcorlor(content)
    data["data"] = content
    data['totalCount'] = totalCount
    
    return JsonResponse(data)

def fetchProteins(request):
    body_unicode = request.body.decode('utf-8')
    body = json.loads(body_unicode)
    print(body)
   
    pageSize = body["pageSize"]
    currentPage = body["currentPage"]
    order = body["order"]
    field = body["field"]
    q_id = body['q_id']

    info = RepeatDomain.objects.get(pk=q_id)
    proteins = RepeatProtein.objects.filter(reprProtein= info.reprProtein)
    
    totalCount = proteins.count() 
    requestData = proteins[(currentPage-1) * pageSize : currentPage * pageSize]
    content = []
    for q in requestData:
        protin_id = q.protein
        cdd_locs, cdd_ids, cdd_notes, cdd_names = [], [], [], []
        protin = NrInfo.objects.filter(protin_id=protin_id).first()
        regions = protCD.objects.filter(protin_id=protin_id)
        for qr in regions:
            qcd = CDD.objects.get(cdd_id=qr.cdd_id_id)
            cdd_ids.append(qr.cdd_id_id)
            cdd_notes.append(qcd.cdd_desc)
            cdd_names.append(qcd.cdd_name)
            cdd_locs.append(
                [qr.start, qr.end, qcd.cdd_name, qcd.cdd_desc])
        # cdd_notes = '^^'.join(cdd_notes)
        items = {
            "q_id":q.id,
            "protin_id": protin_id,
            "length": protin.length,
            "cdd_names": cdd_names,
            "cdd_notes": cdd_notes,
            "cdd_ids": cdd_ids,
            "cdd_locs": cdd_locs,
           
        }
        content.append(items)
    
    data = {"status": 200}
    print("addCDcorlor")
    content = myFunctions.addCDcorlor(content)
    data["data"] = content
    data['totalCount'] = totalCount
    
    return JsonResponse(data)


#
def slim_domain(domains):
    new = []
    for dm in domains:
        if dm['value']:
            new.append(dm)
    return new

def filter_archi(body, data):
    dt = []
    print(body)
    for mydt in data:
        checkDomain = [True]
        for domain in body['domains']:
            if domain['value']:
                isInclude = False
                dCount = domain['count']
                keyword = domain['value'].strip()
                if not domain["exclude"]:
                    # print("key word", keyword, "count: ", dCount,mydt[ domain['type']] )
                    ire_key = f"({keyword}.*?)" + "{" + str(dCount) + ",}"
                    re_k = re.compile(ire_key, re.I)
                    match = re_k.search(mydt[domain['type']])
                    if match:
                        # print('match',match.groups())
                        isInclude = True
                else:
                    # print('aa')
                    re_exclude = re.compile(keyword, re.I)
                    if not re_exclude.search(mydt[domain['type']]):
                        isInclude = True
                checkDomain.append(isInclude)
        # print(checkDomain)
        if not False in checkDomain:
            dt.append(mydt)
    return dt