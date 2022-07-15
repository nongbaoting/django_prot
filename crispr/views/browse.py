from django.http import JsonResponse, FileResponse, Http404

from crispr.models import CASInfo, AlignSPScore, AlignTMScore
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from os import path
import os


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
    objs = AlignTMScore.objects.all()
    totalCount = objs.count()
    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    requestData = objs[(currentPage-1) * pageSize: currentPage * pageSize]
    data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data}
    return JsonResponse(content)


def alignSPscore(request):
    objs = AlignSPScore.objects.all()
    totalCount = objs.count()
    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    requestData = objs[(currentPage-1) * pageSize: currentPage * pageSize]
    data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data}
    return JsonResponse(content)
