pdb_annotations_dir = "/dat1/nbt2/proj/21-prot/web/data/res/pdb_annotations/"
import uuid,time,re,os,json
from protein.toolkit import *
from protein.models import *
from protein import tasks
reg_zip = re.compile('zip$')
reg_W = re.compile("\s+")
from django.http import JsonResponse,FileResponse, Http404
def uploadPDB_and_annotation(request):
    reg_af = re.compile('AF-')
    if request.method == "POST":
        
        fileDict = request.FILES.items()
        # 获取上传的文件，如果没有文件，则默认为None
        if not fileDict:
            return JsonResponse({'msg': 'No file upload'})
        the_date = time.strftime("%Y-%m-%d", time.localtime()) 
        myuuid = "file_" + the_date +'_'+ str(uuid.uuid4())
        tempDir = os.path.join(pdb_annotations_dir, 'uploads',myuuid)
        os.system('mkdir -p ' + tempDir)
        pdbfile= ''
        for (k, v) in fileDict:
            print("dic[%s]=%s" % (k, v))
            fileData = request.FILES.getlist(k)
            for file in fileData:
                fileName = file._get_name()

                fileName = reg_W.sub('_', fileName)
                # filePath = os.path.join(tempDir, fileName)
                filePath = os.path.join(tempDir, 'upload.pdb')
                print('filepath = [%s]' % filePath)
                try:
                    writeFile(filePath, file)
                    pdbfile = filePath
                except:
                    return JsonResponse({'msg': 'file write failed'})

                if reg_zip.search(fileName):
                    myzipfile = zipfile.ZipFile(filePath)
                    # 解压
                    myzipfile.extractall(tempDir)
        
        # params = myFunctions.load_POST(request)
        params ={}
        params['job_name'] = request.POST.get('job_name')
        params['chain'] = request.POST.get('chain')
        print(tempDir)
        params['proj_type'] = "PDB Domain Annotations"
        params['work_dir'] = os.path.join(pdb_annotations_dir,'results')
        params['pdbfile'] =  filePath
        print(params)
        res = tasks.pdb_domain_annotations.apply_async(args = [params])
        myFunctions.create_submit_form(params, res)
        return JsonResponse({'msg':200})


def parser_results(request):
    result_dir = os.path.join(pdb_annotations_dir, 'results')
    myuuid = request.GET.get('uuid')
    request_type  = request.GET.get('request_type')
    if request_type == "interpro":
        resFi = os.path.join(result_dir, myuuid, "interpro.track.json")
        f = open(resFi)
        jsonRes = json.loads(json.load(f))
        return JsonResponse(jsonRes['upload'],safe=False)
    elif request_type == "unidoc":
        resFi = os.path.join(result_dir, myuuid, "unidoc.track.json")
        f = open(resFi)
        jsonRes = json.loads(json.load(f))
        return JsonResponse(jsonRes['upload'],safe=False)
    elif request_type == "sword2":
        resFi = os.path.join(result_dir, myuuid, "sword2.track.json")
        print(request_type, resFi)
        f = open(resFi)
        jsonRes = json.loads(json.load(f))
        return JsonResponse(jsonRes['upload'],safe=False)
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

import numpy as np
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


################################
def writeFile(filePath, file):
    with open(filePath, "wb") as f:
        if file.multiple_chunks():
            for content in file.chunks():
                f.write(content)
        else:
            data = file.read()  # .decode('utf-8')
            f.write(data)