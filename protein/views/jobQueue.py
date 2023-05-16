# from protein.models import *
from collections import defaultdict
import django,json,re
from django.http import JsonResponse
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

from django_celery_results.models import TaskResult
from protein.models import SubmitInfoNew

def jobTable(request):
    pageSize = int(request.GET.get("pageSize"))
    currentPage = int(request.GET.get("currentPage"))
    field = 'job_name' 
    if "field" in request.GET:
        field = request.GET.get("field")
        print("field:", field)
        if request.GET.get("order") == "descending":
            field= '-'+ field
    # 
    taskinfo = TaskResult.objects.all()
    info = SubmitInfoNew.objects.all().order_by('-upload_date')
    totalCount = info.count() 
    # # 分页
    # paginator = Paginator(info, pageSize)
    # print("currentPage, pageSize",currentPage, pageSize)
    # try:
    #     books = paginator.page(currentPage)
    # except PageNotAnInteger:
    #     books = paginator.page(1)
    # except EmptyPage:
    #     books = paginator.page(paginator.num_pages)
    # data  = serializers.serialize('json',books)
    re_com = re.compile(r"'")
    requestData = info[(currentPage-1) * pageSize : currentPage * pageSize]
    newdata=[]
    for subItem in requestData:
        q  = taskinfo.filter(task_id = subItem.uuid).first()
        if not q : continue
        subItem.task_status = q.status
        if not subItem.running_date or not subItem.completed_date:
            subItem.running_date= q.date_created
            subItem.task_status = q.status
            if q.status == "SUCCESS":
                subItem.completed_date= q.date_done
                subItem.result= q.result
        subItem.save()

    data = serializers.serialize('json', requestData)
    data = {"totalCount": totalCount,
            "data": data} 
    return JsonResponse(data)