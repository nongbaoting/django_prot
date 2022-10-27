import json
from collections import defaultdict
from django.http import JsonResponse, FileResponse, Http404
from crispr.models import *
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from os import path
import os
import re
from django.forms.models import model_to_dict

candidateFi = "/dat1/nbt2/proj/22-cas/work/cas9/findCasCrispr/top100/candidate.txt"
candidates = [i for i in open(candidateFi).read().strip(
    '\n').split('\n') if not re.match("#", i)]


def browse(request):
    data = {"hello": "world"}
    obj = CASInfo.objects.all()
    totalCount = obj.count()
    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))

    requestData = obj[(currentPage-1) * pageSize: currentPage * pageSize]
    data = []
    for item in requestData:
        Fi = path.join(
            '/training/nong/protein/work/cas9_ColabFold/colabFoldPdb', f'{item.accession}-cf.pdb')
        pdbExist = path.exists(Fi)
        dt = {
            "cas_class": item.cas_class,
            "accession": item.accession,
            "entry_name": item.entry_name,
            "data_class": item.data_class,
            "protein_name": item.protein_name,
            "gene_name": item.gene_name,
            "organism": item.organism,
            "taxonomy_id": item.taxonomy_id,
            "sequence_length": item.sequence_length,
            "pdb": pdbExist,
        }
        data.append(dt)

    # data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data
               }
    return JsonResponse(content)


re_field = re.compile(r'fields.')


def struc_getFile(request):
    struc_result_dir = "/training/nong/protein/work/cas9_ColabFold/colabFoldPdb"
    if request.method == "GET":
        filename = request.GET.get('filename') + '-cf.pdb'
        filePath = os.path.join(struc_result_dir, filename)
        print(filePath)
        if os.path.exists(filePath):
            print("file ok: ", filePath)
            fo = open(filePath, 'rb')
            return FileResponse(fo)
        else:
            print("not found: ", filePath)
            raise Http404

cas9Dt = {
    'spCas9-3': ['4un3'],
    'cjCas9-3': ['5x2g',"AF-Q0P897-F1-model_v3"],
    'Nme1Cas9-3': ['6jdv']
}

def alignTMscore(request):

    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    field = '-seq_ID'
    if request.GET.get("field"):
        # print(request.GET.get('field'))
        field = re_field.sub('', request.GET.get("field"))
        # print("field:", field)
        if request.GET.get("order") == "descending":
            field = '-' + field

    filters = json.loads(request.GET.get('filters'))
    objs = AlignTMScore.objects.all().filter(
        chain1__in=cas9Dt[filters['protein']],
        chain2_len__gt=filters['min_len'], chain2_len__lt=filters['max_len'],
        seq_ID__gt=filters['min_SI'], seq_ID__lt=filters['max_SI']
    ).order_by(field)

    if filters['exclude_knownCas'] == True:
        objs = filter_known_cas(objs)
    if filters['candidates'] == True:
        objs = objs.filter(chain2_acc__in=candidates)
    totalCount = objs.count()
    requestData = objs[(currentPage-1) * pageSize: currentPage * pageSize]

    data = []
    for item in requestData:
        it = model_to_dict(item)
        # print(it)
        cas9 = CASInfo.objects.get(accession=item.chain2_acc)
        it["protein_name"] = cas9.protein_name
        it["organism"] = cas9.organism
        it["genome_genbank"] = cas9.genome_genbank
        it["protein_genebankID"] = cas9.protein_genebankID
        data.append(it)
    # data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data}
    return JsonResponse(content)


def alignSPscore(request):
    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    tool = request.GET.get('tool')
    field = '-seq_ID'
    filters = json.loads(request.GET.get('filters'))
    # print(filters)
    if request.GET.get("field"):
        print(request.GET.get('field'))
        field = re_field.sub('', request.GET.get("field"))
        print("field:", field)
        if request.GET.get("order") == "descending":
            field = '-' + field
    objs = AlignSPScore.objects.all().filter(
        chain1__in=cas9Dt[filters['protein']],
        tool=tool).all().filter(chain2_len__gt=filters['min_len'], chain2_len__lt=filters['max_len'], seq_ID__gt=filters['min_SI'], seq_ID__lt=filters['max_SI']).order_by(field)

    if filters['exclude_knownCas'] == True:
        objs = filter_known_cas(objs)
    if filters['candidates'] == True:
        objs = objs.filter(chain2_acc__in=candidates)

    totalCount = objs.count()
    print('-----------', totalCount)
    requestData = objs[(currentPage-1) * pageSize: currentPage * pageSize]
    data = []
    for item in requestData:
        it = model_to_dict(item)
        # print(it)
        cas9 = CASInfo.objects.get(accession=item.chain2_acc)
        it["protein_name"] = cas9.protein_name
        it["organism"] = cas9.organism
        it["genome_genbank"] = cas9.genome_genbank
        it["protein_genebankID"] = cas9.protein_genebankID
        data.append(it)
    # data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data}
    return JsonResponse(content)


def alignFatcatScore(request):
    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))

    tool = request.GET.get('tool')
    field = '-seq_ID'
    filters = json.loads(request.GET.get('filters'))
    print(filters)
    if request.GET.get("field"):
        print(request.GET.get('field'))
        field = re_field.sub('', request.GET.get("field"))
        print("field:", field)
        if request.GET.get("order") == "descending":
            field = '-' + field

    objs = AlignFatcatScore.objects.all().filter(
        chain1__in=cas9Dt[filters['protein']],).all().filter(chain2_len__gt=filters['min_len'],
                                                         seq_ID__gt=filters['min_SI'], seq_ID__lt=filters['max_SI'],
                                                         chain2_len__lt=filters['max_len']).order_by(field)

    if filters['exclude_knownCas'] == True:
        objs = filter_known_cas(objs)

    if filters['candidates'] == True:
        objs = objs.filter(chain2_acc__in=candidates)
    totalCount = objs.count()
    print('-----------', totalCount)
    requestData = objs[(currentPage-1) * pageSize: currentPage * pageSize]
    data = []
    for item in requestData:
        it = model_to_dict(item)
        # print(it)
        cas9 = CASInfo.objects.get(accession=item.chain2_acc)
        it["protein_name"] = cas9.protein_name
        it["organism"] = cas9.organism
        it["genome_genbank"] = cas9.genome_genbank
        it["protein_genebankID"] = cas9.protein_genebankID
        it['repeatinfo'] = repeatDt[cas9.genome_genbank]
        data.append(it)
    # data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data}
    return JsonResponse(content)


########### functions #########
infoFi = "/dat1/nbt2/proj/22-cas/work/cas9/known_cas9.txt"
info = open(infoFi).read().strip('\n').split('\n')
cas9Info = defaultdict(list)
for li in info:
    cell = li.split("\t")
    cas9Info[cell[1]] = cell

def filter_known_cas(scoreObj):
    casinfo = CASInfo.objects.all()
    for cas9 in cas9Info:
        casinfo = casinfo.exclude(organism__contains=cas9Info[cas9][2])
    scoreObj = scoreObj.filter(chain2_acc__in=casinfo.values_list(
        'accession', flat=True))
    return scoreObj


repeatFi = "/dat1/nbt2/proj/22-cas/work/cas9/findCasCrispr/top100/results/anti_repeat.info"
repeatDt = defaultdict(list)
with open(repeatFi) as f:
    for li in f:
        cell = li.strip('\n').split("\t")
        repeatDt[cell[0] ].append  ( cell)

