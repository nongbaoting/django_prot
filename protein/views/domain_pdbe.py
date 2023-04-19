from django.http import JsonResponse, HttpResponse, FileResponse, Http404
from protein.toolkit import interproscan
from django.core import serializers
import json
def interpro(request):
    fi = "/dat1/nbt2/proj/23-tadA/work/temp/CDKAL_HUMAN.json"
    inter = interproscan.InterproScan()
    data = inter.parse(fi)
    item = data['CDKAL_HUMAN']
    # rowConfigData = serializers.serialize('json', data)
    
    return JsonResponse(item)