import os, pickle
from os import path
import re
import zipfile
import uuid,json
from pydoc import describe
from django.conf import settings
from django.http import JsonResponse,FileResponse, Http404
from collections import defaultdict
from protein.toolkit import *
from protein.models import *
from protein import tasks
from py4j.java_gateway import JavaGateway
from protein.toolkit import TMalign
# jarpath = "/training/nong/web/java/CalProSimilariry-1.0.1-jar-with-dependencies.jar"
# 读入biozernike 信息
human_biozernike = "/training/nong/protein/biozernike/human/human_biozernike.json"
pdb_biozernike = "/dat1/nbt2/proj/21-prot/dat/pdb/biozernike/pdb_biozernike.json"
gateway = JavaGateway()
mycompare = gateway.entry_point

# human = mycompare.loadDesciptor(human_biozernike)
# pdb = mycompare.loadDesciptor(pdb_biozernike)
# dev save
human = ''
pdb = ''
biozernike = {
    "human": human,
    'pdb': pdb
}
reg_zip = re.compile('zip$')
reg_W = re.compile("\s+")
struc_cpm_dir = "/dat1/nbt2/proj/21-prot/web/data/res/structure_comparison"
scopeDomain_dir ="/dat1/nbt2/proj/21-prot/dat/pdb/scope_domain"

################ view ###############

def aligment(request):
    if request.method == "GET":
        pdb1 = request.GET.get('pdb1')
        pdb2 = request.GET.get('pdb2')
        pdb1_file=''
        pdb2_file=''
        if path.exists(pdb1_file) and path.exists(pdb2_file):
            pass

def DUF_SPalign(request):
    jsonFi = "/training/nong/web/data/res/struc_comparison/duf_q2_cov.6.json"
    with open(jsonFi)  as f:
        data_str = json.load(f)
        data = {"data":data_str}
        return JsonResponse(data)


def results(request):
    if request.method =="GET":
        myuuid = request.GET.get('uuid')
        resFi = os.path.join(struc_cpm_dir, myuuid, 'TMalgin.pickle')
        print(resFi)
        res = myFunctions.pickle_load_file(resFi)
        dt = get_one_page(res, request)
        data = {'data':dt,
        'totalCount':len(res)}
        return JsonResponse(data)

def getOneItem(request):
    if request.method =="GET":
        myuuid = request.GET.get('uuid')
        dataType = request.GET.get("dataType")
        tmpdir = os.path.join(struc_cpm_dir, myuuid)
        input_pdb = os.path.join(tmpdir, request.GET.get('input_pdb') )
        db_pdb =  os.path.join(scopeDomain_dir, request.GET.get('db_pdb') )
        outFi_name = request.GET.get('input_pdb') + request.GET.get('db_pdb')
        item_pickle = os.path.join(tmpdir, outFi_name + '.pickle')
        if dataType == 'info':
            if not os.path.exists(item_pickle):
                item = TMalign.one_to_one(input_pdb, db_pdb, tmpdir,outFi_name)
            else:
                item = myFunctions.pickle_load_file(item_pickle)
            return JsonResponse(add_scopecla(item))
        if dataType == 'input_pdb':
            pdb_fi = os.path.join(tmpdir, outFi_name + '.sup.pdb')
            print(pdb_fi)
            fh = open(pdb_fi, 'rb')
            return FileResponse(fh)
        elif dataType == "db_pdb":
            pdb_fi = db_pdb
            if os.path.exists(pdb_fi):
                print(pdb_fi)
                fh = open(pdb_fi, 'rb')
                return FileResponse(fh)
            else:
                print("File not exist!" + pdb_fi)
                raise Http404("File not exist!" + pdb_fi)

def upload_pdb(request):
    reg_af = re.compile('AF-')
    if request.method == "POST":
        fileDict = request.FILES.items()
        # 获取上传的文件，如果没有文件，则默认为None
        if not fileDict:
            return JsonResponse({'msg': 'No file upload'})

        tempDirName = "file_" + str(uuid.uuid4())
        tempDir = os.path.join(struc_cpm_dir, 'tmp', tempDirName)
        os.system('mkdir -p ' + tempDir)

        for (k, v) in fileDict:
            print("dic[%s]=%s" % (k, v))
            fileData = request.FILES.getlist(k)
            for file in fileData:
                fileName = file._get_name()

                fileName = reg_W.sub('_', fileName)
                filePath = os.path.join(tempDir, fileName)
                print('filepath = [%s]' % filePath)
                try:
                    writeFile(filePath, file)
                except:
                    return JsonResponse({'msg': 'file write failed'})

                if reg_zip.search(fileName):
                    myzipfile = zipfile.ZipFile(filePath)
                    # 解压
                    myzipfile.extractall(tempDir)
        
        # params = myFunctions.load_POST(request)
        params ={}
        params['dir_1'] = tempDir
        params['dir_2'] = scopeDomain_dir
        params['job_name'] =request.POST.get('job_name')
        print(tempDir)
        params['proj_type'] = "Structure Comparison"
        params['struc_cpm_dir'] = struc_cpm_dir
        params['outFi_name'] = "TMalgin.pickle"
        res = tasks.structure_comparison.apply_async(args = [params])
        myFunctions.create_submit_form(params, res)
        return JsonResponse({'msg':200})

################## function ###############################

def get_one_page(newQ, request):
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    order="descending"
    field= "tmscore_1"
    order = request.GET.get('order')
    isSort = request.GET.get('isSort')
    field = request.GET.get('field')
    if isSort:
        print("order field", order, field,isSort)
        if order == "descending":
            newQ = sorted(newQ, key=lambda k: k[field], reverse=True)
        else:
            newQ = sorted(newQ, key=lambda k: k[field])

    requestData = newQ[(currentPage - 1) * pageSize : currentPage * pageSize]
    content = []
    for item in requestData:
        add_scopecla(item)
        content.append( add_scopecla(item) )
    return  content


def add_scopecla(item):
    domain_id = item['chain_2_scopeDomain']
    obj = ScopeCla.objects.filter(domain = domain_id).first()
    item['family'] = obj.family
    item['superfamily'] = obj.superfamily
    return item

def scanAndFind_notZip(mydir):
    wantdir = []
    for entry in os.scandir(mydir):
        if entry.is_file() and not reg_zip.search(entry.name):
            wantdir.append(entry)
        elif entry.is_dir():
            wantdir.extend(scanAndFind_notZip(entry.path))
    return wantdir

def writeFile(filePath, file):
    with open(filePath, "wb") as f:
        if file.multiple_chunks():
            for content in file.chunks():
                f.write(content)
        else:
            data = file.read()  # .decode('utf-8')
            f.write(data)


def compare_descriptor(pdbfile, species):
    # pdb = "/training/nong/protein/biozernike/human/one_test/test/1e00.cif.gz"
    desc = biozernike[species]
    scores = mycompare.compare(pdbfile, desc)
    scores_cut = mycompare.sortByValueArray_cut(scores)
    # print(scores_cut)
    scDict = dict(scores_cut)
    return scDict
