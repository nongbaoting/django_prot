pdb_annotations_dir = "/dat1/nbt2/proj/21-prot/web/data/res/pdb_annotations/"
import uuid,time,re,os,json
from protein.toolkit import *
from protein.models import *
from protein import tasks
reg_zip = re.compile('zip$')
reg_W = re.compile("\s+")
reg_cif = re.compile('cif$')

from Bio import PDB
import numpy as np

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
                if reg_cif.search(fileName):
                    filePath = os.path.join(tempDir, 'upload.cif')
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
        with open(resFi) as f:
            jsonRes = json.loads(json.load(f))
        return JsonResponse(jsonRes['upload:A'],safe=False)
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
            # jsonRes = json.loads(json.load(f))
        # print(jsonRes)
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



def align(request):
    myuuid = request.GET.get('uuid')
    db_pdbid =request.GET.get('db_pdbid')
    db_name = request.GET.get('db_name')
    chain = request.GET.get('chain')
    resFi = os.path.join(pdb_annotations_dir, 'results', myuuid, "ecod.json")
    print(request)
    db_dir = '/dat1/dat/db/ECOD/F70/data/ecod/domain_data_ln/'
    suffix = '.pdbnum.pdb'
    db_pdb = os.path.join(db_dir,db_pdbid + suffix)
    upload_pdb = os.path.join(pdb_annotations_dir, 'results', myuuid, "upload.pdb")
    with open(resFi) as f:
        jsonRes = json.loads(json.load(f))
        matrix = getMatrix(db_pdbid, jsonRes)
        print(matrix)
        newUploadPdb = rotate(matrix, upload_pdb)
    rotate_pdb = "/dat1/nbt2/proj/21-prot/web/data/res/pdb_annotations/results/e1960a66-ea52-48f0-905b-f4a4f6fc201f/test.pdb"
    
    parser = PDB.PDBParser()
    struct1 = parser.get_structure('protein1', rotate_pdb)

    # Read second PDB file
    struct2 = parser.get_structure('protein2', db_pdb)
    for model in struct2:
        for chain in model:
            chain.id = 'B' # Rename chain ID
            struct1[0].add(chain) 
    # Merge the two structures
    # merged = struct1 + struct2

    # Write out merged structure to StringIO buffer
    # io = BytesIO()
    # PDB.PDBIO().set_structure(struct1).save(io)
    # io.seek(0)

    # response = FileResponse(struct2, filename='merged.pdb', content_type='chemical/x-pdb')
    # return response
    fh = open(rotate_pdb, 'rb')
    return FileResponse(fh)



################################

def getMatrix(target, alignRes):
    matrix = []
    for item in alignRes:
        if item["target"] == target:
            t = [float(i) for i in item['t'].split(',')]
            u = [float(i) for i in item['u'].split(',')]
            print(t)
            matrix =  np.array([
                [t[0], u[0], u[3],u[6]],
                [t[1], u[1], u[4],u[7]],
                [t[2], u[2], u[5],u[8]],
                ])
            return matrix

def rotate(matrix, inPDB_1):
        parser = PDB.PDBParser()
        name = inPDB_1.split('.')[0]
        struct = parser.get_structure(name, inPDB_1)
        for model in struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        x, y, z = atom.coord
                        X = matrix[0][0]     + matrix[0][1] * x + \
                            matrix[0][2] * y + matrix[0][3] * z

                        Y = matrix[1][0]     + matrix[1][1] * x + \
                            matrix[1][2] * y + matrix[1][3] * z

                        Z = matrix[2][0]     + matrix[2][1] * x + \
                            matrix[2][2] * y + matrix[2][3] * z

                        atom.coord = [X, Y, Z]
        io = PDB.PDBIO()
        io.set_structure(struct)
        io.save('/dat1/nbt2/proj/21-prot/web/data/res/pdb_annotations/results/e1960a66-ea52-48f0-905b-f4a4f6fc201f/test.pdb')
        return struct


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