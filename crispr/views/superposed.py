import uuid
from django.http import JsonResponse, FileResponse, Http404

from crispr.models import CASInfo, AlignSPScore, AlignTMScore
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from os import path
import os
import re
from crispr.utils.SPalign import run_spalign

process_dir = "/dat1/nbt2/proj/21-prot/web/data/res/structure_comparison"


def getAlignPDBFile(request):
    tool = request.GET.get('tool')
    filename = request.GET.get('filename')
    uuid = request.GET.get('uuid')
    filePath = ''
    if tool == 'TMalign':
        filePath = os.path.join(filePath, uuid, tool, filename)
    fh = open(filePath, 'rb')
    return FileResponse(fh, filename="aa")


def alignment(request):
    tool = request.GET.get('tool')
    input_pdb = request.GET.get('input_pdb')
    target_pdb = request.GET.get('target_pdb')
    content = {}
    tempDir, uuid_dir = create_tmpDir(process_dir, request)
    if tool == 'TMalign':
        pass
    if tool == 'SPalign':
        print("Spalgin")
        spalign = run_spalign(input_pdb, target_pdb, tempDir)
        content = spalign.content
    return JsonResponse(content)

# logic


def create_tmpDir(process_dir, info):
    tool = info.GET.get('tool')
    input_pdb = path.basename(info.GET.get('input_pdb'))
    target_pdb = path.basename(info.GET.get('target_pdb'))
    uuid_dir = str(uuid.uuid1())
    tempDir = '.'.join([uuid_dir, tool, input_pdb, target_pdb])
    tempDir = path.join(process_dir, tempDir)
    os.system("mkdir -p " + tempDir)
    return [tempDir, uuid_dir]
