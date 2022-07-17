from django.http import JsonResponse, FileResponse, Http404

from crispr.models import CASInfo, AlignSPScore, AlignTMScore
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from os import path
import os
import re


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
    pass
