import netifaces 
pdb_annotations_dir = "/apps/resultData/pdb_annotations/"
pdb_annotations_dir_remote = "/home/public/nong/pdb_annotations/"

ip=''
import uuid,time,re,os,json
from protein.toolkit import *
from protein.models import *
import prody
from protein import tasks
reg_zip = re.compile('zip$')
reg_W = re.compile("\s+")
reg_cif = re.compile('cif$')
from Bio import PDB
import numpy as np
from django.http import JsonResponse,FileResponse, Http404
from django.core import serializers
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

def browse(request):
    request_type = request.GET.get('request_type')
    if request_type == "browse":
        pageSize = int(request.GET.get("pageSize"))
        currentPage = int(request.GET.get("currentPage"))
        #uniprot_id = request.GET.get('uniprot_id')
        info = AF_uniprot.objects.all()
        totalCount = info.count() 
        # 分页
        paginator = Paginator(info, pageSize)
        print("currentPage, pageSize",currentPage, pageSize)
        try:
            dataPage = paginator.page(currentPage)
        except PageNotAnInteger:
            dataPage = paginator.page(1)
        except EmptyPage:
            dataPage = paginator.page(paginator.num_pages)
        data  = serializers.serialize('json',dataPage)
        return JsonResponse({"data":data, "totalCount":totalCount }, safe=False)
    elif request_type == "search":
        pageSize = int(request.GET.get("pageSize"))
        currentPage = int(request.GET.get("currentPage"))
        #uniprot_id = request.GET.get('uniprot_id')
        searchType = request.GET.get("searchType")
        searchContent = request.GET.get('searchContent')
        filters = searchType + '__icontains'
        print(searchType, searchContent,filters)
        if searchType=="gene_names":
            info = AF_uniprot.objects.filter(gene_names__icontains = searchContent)
        elif searchType=="protein_names":
            info = AF_uniprot.objects.filter(protein_names__icontains = searchContent)
        elif not searchContent:
            info = AF_uniprot.objects.all()
        # info = AF_uniprot.objects.filter(**{filters: searchContent})
        totalCount = info.count()
        # 分页
        paginator = Paginator(info, pageSize)
        print("currentPage, pageSize",currentPage, pageSize)
        try:
            dataPage = paginator.page(currentPage)
        except PageNotAnInteger:
            dataPage = paginator.page(1)
        except EmptyPage:
            dataPage = paginator.page(paginator.num_pages)
        data  = serializers.serialize('json',dataPage)
        # data  = serializers.serialize('json',info)
        return JsonResponse({"data":data, "totalCount":totalCount }, safe=False)
    
def uploadPDB_and_annotation(request):
    reg_af = re.compile('AF-')
    if request.method == "POST":
        fileDict = request.FILES.items()
        # 获取上传的文件，如果没有文件，则默认为None
        if not fileDict:
            return JsonResponse({'msg': 'No file upload'})
        the_date = time.strftime("%Y-%m-%d", time.localtime()) 
        myuuid = "file_" + the_date +'_'+ str(uuid.uuid4())
        uploadDir = os.path.join(pdb_annotations_dir, 'uploads',myuuid)
        os.system('mkdir -p ' + uploadDir)
        pdbfile= ''
        for (k, v) in fileDict:
            print("dic[%s]=%s" % (k, v))
            fileData = request.FILES.getlist(k)
            for file in fileData:
                fileName = file._get_name()

                fileName = reg_W.sub('_', fileName)
                # filePath = os.path.join(uploadDir, fileName)
                # TODO 处理用户上传 chain 选择
                filePath = os.path.join(uploadDir, 'uploadRaw.pdb')
                if reg_cif.search(fileName):
                    filePath = os.path.join(uploadDir, 'uploadRaw.cif')
                print('filepath = [%s]' % filePath)
                try:
                    writeFile(filePath, file)
                    pdbfile = filePath
                except:
                    return JsonResponse({'msg': 'file write failed'})

                if reg_zip.search(fileName):
                    myzipfile = zipfile.ZipFile(filePath)
                    # 解压
                    myzipfile.extractall(uploadDir)
        
        
        params ={}
        params['job_name'] = request.POST.get('job_name')
        params['chain'] = request.POST.get('chain')
        print(uploadDir)
        params['proj_type'] = "PDB Domain Annotations"
        params['work_dir'] = os.path.join(pdb_annotations_dir, 'results')
        params['pdb_annotations_dir_remote'] = pdb_annotations_dir_remote
        params['pdbfile'] =  pdbfile
        params['ip'] = ip
        print(params)
        res = tasks.pdb_domain_annotations.apply_async(args = [params])
        
        myFunctions.create_submit_form(params, res)
        content = {'msg': 200,
                   'uuid': res.task_id
                   
                   }
        print(content)
        return JsonResponse(content)

def check_job(request):
    uuid = request.GET.get('uuid')
    content = {
        'uuid': uuid,
        'task_status': 'not_exist'
    }
    job = SubmitInfoNew.objects.filter(uuid=uuid).first()
    if job:
        content['task_status'] = job.task_status
    return JsonResponse(content)

def parser_results(request):
    result_dir = os.path.join(pdb_annotations_dir, 'results')
    myuuid = request.GET.get('uuid')
    request_type  = request.GET.get('request_type')
    resFi= ''
    if request_type == "datainfo":
        Fi = os.path.join(result_dir, myuuid, "params.json")
        with open(Fi) as f:
            params = json.loads(json.load(f))
            content = {'msg':200, 'chain':params['chain']}
        return  JsonResponse(content)
    elif request_type == "interpro":
        resFi = os.path.join(result_dir, myuuid, "interpro.track.json")
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            return JsonResponse(jsonRes['upload'],safe=False)
    elif request_type == "unidoc":
        resFi = os.path.join(result_dir, myuuid, "unidoc.track.json")
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            return JsonResponse(jsonRes['upload'],safe=False)
    elif request_type == "sword2":
        resFi = os.path.join(result_dir, myuuid, "sword2.track.json")
        print(request_type, resFi)
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            return JsonResponse(jsonRes['upload'],safe=False)
    elif request_type == "protvista":
        resFi = os.path.join(result_dir, myuuid, "protvistData.json")
        print(request_type, resFi)
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            return JsonResponse(jsonRes,safe=False)
    elif request_type == "ECOD":
        resFi = os.path.join(result_dir, myuuid, "ecod.json")
        print(request_type, resFi)
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            data = get_one_page(jsonRes,request)
            content = {
                "data": data,
                "totalCount": len(jsonRes)
            }
            return JsonResponse(content,safe=False)
    elif request_type == "SCOP":
        resFi = os.path.join(result_dir, myuuid, "scop.json")
        print(request_type, resFi)
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            data = get_one_page(jsonRes,request)
            content = {
                "data": data,
                "totalCount": len(jsonRes)
            }
            return JsonResponse(content,safe=False)
    elif request_type == "CATH":
        resFi = os.path.join(result_dir, myuuid, "cath.json")
        print(request_type, resFi)
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            data = get_one_page(jsonRes,request)
            content = {
                "data": data,
                "totalCount": len(jsonRes)
            }
            return JsonResponse(content,safe=False)
    elif request_type == "AFDB":
        resFi = os.path.join(result_dir, myuuid, "AFDB.json")
        print(request_type, resFi)
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            data = get_one_page(jsonRes,request)
            content = {
                "data": data,
                "totalCount": len(jsonRes)
            }
            return JsonResponse(content,safe=False)
    elif request_type == "pdbDB":
        resFi = os.path.join(result_dir, myuuid, "pdbDB.json")
        print(request_type, resFi)
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
            data = get_one_page(jsonRes,request)
            content = {
                "data": data,
                "totalCount": len(jsonRes)
            }
            return JsonResponse(content,safe=False)
    else:
        raise Exception("Unknown")

    
def get_pdbFile(request):
    result_dir = os.path.join(pdb_annotations_dir, 'results')
    myuuid = request.GET.get('uuid')
    pdb_fi = os.path.join(result_dir, myuuid, "upload.pdb")
    if os.path.exists(pdb_fi):
        print(pdb_fi)
        fh = open(pdb_fi, 'rb')
        return FileResponse(fh)
    else:
        print("File not exist!" + pdb_fi)
        raise Http404("File not exist!" + pdb_fi)
def get_alignPDBfile(request):
    pdb_fi = request.GET.get("pdbFile")
    if os.path.exists(pdb_fi):
        print(pdb_fi)
        fh = open(pdb_fi, 'rb')
        return FileResponse(fh)
    else:
        print("File not exist!" + pdb_fi)
        raise Http404("File not exist!" + pdb_fi)

def align(request):
    annotateDBinfo = {
        "ECOD": ['/apps/data/PDB/ECOD/'],
        "SCOP": ['/apps/data/PDB/SCOP/'],
        "pdbDB": ["/apps/data/PDB/pdbDB"],
        "AFDB": ["/apps/data/PDB/AFDB"],
        "CATH":[]
    }
    
    myuuid = request.GET.get('uuid')
    db_pdbName =request.GET.get('db_pdbid').split('_MODEL')[0]
    db_name = request.GET.get('db_name')
    chain = request.GET.get('chain')
    work_dir = os.path.join(pdb_annotations_dir, 'results', myuuid)
    # resFi = os.path.join(work_dir, "ecod.json")
    print(request)
    db_dir = annotateDBinfo[db_name][0]
    db_pdb = os.path.join(db_dir, db_pdbName)
    print("db_pdb:", db_pdb)
    if db_name in ['pdbDB']:
        db_pdb = process_dbPDB(work_dir, db_pdb)
    
    upload_pdb = os.path.join(work_dir, "upload.pdb")
    # TODO
    # tmalign = TMalign(work_dir, upload_pdb, db_pdb)
    # tmalign.run()
    params = {}
    params['work_dir'] = work_dir
    params['upload_pdb'] = upload_pdb
    params['db_pdb'] = db_pdb
    
    res = tasks.Alignment.apply_async(args = [params] )
    
    content = {'task_id':res.id}
    return JsonResponse(content)
    # print(TMAlign_mergePDB)
    # fh = open(TMAlign_mergePDB, 'rb')
    # return FileResponse(fh)

from celery.result import AsyncResult
def check_TMalign(request):
    task_id = request.GET.get('task_id')
    async_result = AsyncResult(task_id)
    try:
        result = async_result.get(timeout=5, propagate=False)
    except TimeoutError:
        result = None
    status = async_result.status
    traceback = async_result.traceback
    if isinstance(result, Exception):
        return JsonResponse({
            'status': status,
            'error': str(result),
            'traceback': traceback,
        })
    else:
        return JsonResponse({
            'status': status,
            'result': result,
        })

def get_alignPDBfile(request):
    pdbFile = request.GET.get('pdbFile')
    if os.path.exists(pdbFile):
        print(pdbFile)
        fh = open(pdbFile, 'rb')
        return FileResponse(fh)
    else:
        print("File not exist!" + pdbFile)
        raise Http404("File not exist!" + pdbFile)
    
################################
re_sp = re.compile(r'_')
def process_dbPDB(work_dir, db_pdb):
    if not re_sp.search(db_pdb):
        return db_pdb
    db_dir = os.path.dirname(db_pdb)
    name = os.path.basename(db_pdb)
    source_pdb_name, chain = re_sp.split(name,1)
    source_pdb = os.path.join(db_dir, source_pdb_name)
    pdbid = source_pdb_name.split('.')[0] + '.pdb'
    new_pdb = os.path.join(work_dir,'tmpAlign',pdbid)
    os.system("mkdir -p " + os.path.join(work_dir,'tmpAlign'))
    atoms = prody.parseMMCIF(source_pdb, chain=chain)
    if atoms:
        writePrody = prody.writePDB(new_pdb, atoms)
    return new_pdb

def example_pdb(request):
    pdbFile = "./protein/static/upload.pdb"
    if os.path.exists(pdbFile):
        print(pdbFile)
        fh = open(pdbFile, 'rb')
        return FileResponse(fh)
    else:
        print("File not exist!" + pdbFile)
        raise Http404("File not exist!" + pdbFile)


def df_to_JSjson(df):
    columns = df.columns
    print(columns)
    data = []
    for row_index,row in df.iterrows():
        item = {}
        for idx,colname in enumerate(list(columns)):
            value = row[idx]
            if not value or value is np.nan:
                value = 0
            item[colname] = value
        data.append(item)
    return data

def tableFi_toJSjson(fi):
    lines = open(fi).read().strip('\n').split('\n')
    headers = lines[0].split('\t')
    dt = []
    for li in lines[1:]:
        cell = li.split('\t')
        item = {}
        for idx, headname in enumerate(headers):
            value =cell[idx]
            if cell[idx] == 'NA':
                value = ''
            item[headname] = value
        dt.append(item)
    return dt


def writeFile(filePath, file):
    with open(filePath, "wb") as f:
        if file.multiple_chunks():
            for content in file.chunks():
                f.write(content)
        else:
            data = file.read()  # .decode('utf-8')
            f.write(data)

def get_one_page(newQ, request, field="ttmscore",order="descending"):
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    
    order = request.GET.get('order')
    is_sort = request.GET.get('is_sort')
    field = request.GET.get('field')
    print("order field", order, field,is_sort)
    if is_sort:
        print("order field", order, field,is_sort)
        if order == "descending":
            newQ = sorted(newQ, key=lambda k: k[field], reverse=True)
        else:
            newQ = sorted(newQ, key=lambda k: k[field])

    requestData = newQ[(currentPage - 1) * pageSize : currentPage * pageSize]
    return  requestData