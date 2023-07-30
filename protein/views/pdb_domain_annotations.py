pdb_annotations_dir = "/dat1/nbt2/proj/21-prot/web/data/res/pdb_annotations/"
import uuid,time,re,os,json
from protein.toolkit import *
from protein.models import *
from protein import tasks
reg_zip = re.compile('zip$')
reg_W = re.compile("\s+")
reg_cif = re.compile('cif$')
import shutil

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
    resFi= ''
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




def align(request):
    annotateDBinfo = {
        "ECOD": ['/dat1/dat/db/ECOD/F70/data/ecod/domain_data_ln/'],
        "SCOP": ['/dat1/nbt2/proj/21-prot/dat/pdb/scope_domain'],
        "pdbDB": ["/dat1/nbt2/proj/21-prot/dat/pdb/pdb_ln"],
        "AFDB": ["/dat1/nbt2/proj/21-prot/dat/alphafold/v2/all/AFDB_all"],
        "CATH":[]
    }   
    myuuid = request.GET.get('uuid')
    db_pdbName =request.GET.get('db_pdbid').split('_MODEL')[0]
    db_name = request.GET.get('db_name')
    chain = request.GET.get('chain')
    work_dir = os.path.join(pdb_annotations_dir, 'results', myuuid)
    resFi = os.path.join(work_dir, "ecod.json")
    print(request)
    db_dir = annotateDBinfo[db_name][0]
    db_pdb = os.path.join(db_dir,db_pdbName)
    upload_pdb = os.path.join(work_dir, "upload.pdb")

    tmalign = TMalign(work_dir, upload_pdb, db_pdb)
    tmalign.run()
    fh = open(tmalign.merged_pdb, 'rb')
    return FileResponse(fh)



################################
class TMalign:
    def __init__(self, work_dir, pdb1, pdb2):
        self.matrixFi = os.path.join(work_dir, 'matrix')
        self.work_dir = work_dir
        self.pdb1 = pdb1 
        self.pdb2 = pdb2
       
        self.merged_pdb = os.path.join(work_dir, 'merged_pdb.pdb')
        self.matrix = []
    def clean(self):
        tmpFiles = [self.matrixFi,self.merged_pdb]
        for tFi in tmpFiles:
            if os.path.exists(tFi):
                os.remove(tFi)
                
    def getMatrix(self):
        re_ma = re.compile(r"m\s+t\[m\]\s+u\[m\]\[0\]")
        re_sp = re.compile(r"\s+")
        matrix = []
        with open(self.matrixFi, 'r') as f:
            m = 0
            for li in f:
                if m == 1:
                    cell = re_sp.split(li.strip())
                    self.matrix.append([float(i) for i in cell[1:]])
                if re_ma.match(li):
                    m += 1
                if len(self.matrix) >= 3:
                    print(self.matrix)
                    break

    def rotate(self):
        parser = PDB.PDBParser()
        self.struct1 = parser.get_structure('uploadPDB', self.pdb1)
        for model in self.struct1:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        x, y, z = atom.coord
                        X = self.matrix[0][0] + self.matrix[0][1] * x + \
                            self.matrix[0][2] * y + self.matrix[0][3] * z

                        Y = self.matrix[1][0] + self.matrix[1][1] * x + \
                            self.matrix[1][2] * y + self.matrix[1][3] * z

                        Z = self.matrix[2][0] + self.matrix[2][1] * x + \
                            self.matrix[2][2] * y + self.matrix[2][3] * z

                        atom.coord = [X, Y, Z]


    def merge2pdb(self,):
        parser = PDB.PDBParser()
        # struct1 = parser.get_structure('protein1', pdb1)
        # Read second PDB file
        struct2 = parser.get_structure('protein2', self.pdb2)
        for model in struct2:
            for chain in model:
                chain.id = 'B' # Rename chain ID
                self.struct1[0].add(chain)
        io = PDB.PDBIO()
        io.set_structure(self.struct1)
        io.save(self.merged_pdb)

    def run_cmd(self):
        cmd = f'cd {self.work_dir}; TMalign {self.pdb1} {self.pdb2} -m {self.matrixFi}'
        print(cmd)
        os.system(cmd)

    def run(self,):
        self.run_cmd()
        self.getMatrix()
        self.rotate()
        self.merge2pdb()

# TODO delete obselete
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
# TODO delete obselete
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