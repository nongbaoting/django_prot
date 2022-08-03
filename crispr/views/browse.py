import json
from django.http import JsonResponse, FileResponse, Http404

from crispr.models import CASInfo, AlignSPScore, AlignTMScore
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from os import path
import os
import re


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
               "data": data}
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


def alignTMscore(request):

    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    field = '-chain2_len'
    if request.GET.get("field"):
        print(request.GET.get('field'))
        field = re_field.sub('', request.GET.get("field"))
        print("field:", field)
        if request.GET.get("order") == "descending":
            field = '-' + field

    filters = json.loads(request.GET.get('filters'))
    objs = AlignTMScore.objects.all().filter(chain2_len__gt=filters['min_len'], chain2_len__lt=filters['max_len'],
                                             seq_ID__gt=filters['min_SI'], seq_ID__lt=filters['max_SI']
                                             ).order_by(field)
    totalCount = objs.count()
    requestData = objs[(currentPage-1) * pageSize: currentPage * pageSize]
    data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data}
    return JsonResponse(content)


def alignSPscore(request):
    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    tool = request.GET.get('tool')
    field = '-chain2_len'
    filters = json.loads(request.GET.get('filters'))
    print(filters)
    if request.GET.get("field"):
        print(request.GET.get('field'))
        field = re_field.sub('', request.GET.get("field"))
        print("field:", field)
        if request.GET.get("order") == "descending":
            field = '-' + field
    objs = AlignSPScore.objects.all().filter(
        tool=tool).all().filter(chain2_len__gt=filters['min_len'], chain2_len__lt=filters['max_len'], SPb__gt=0.5).order_by(field)
    totalCount = objs.count()
    print('-----------', totalCount)
    requestData = objs[(currentPage-1) * pageSize: currentPage * pageSize]
    data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data}
    return JsonResponse(content)
