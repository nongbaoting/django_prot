import os
import re
import zipfile
import uuid,json
from pydoc import describe
from django.conf import settings
from django.http import JsonResponse
from collections import defaultdict
from protein.toolkit import *
from protein import tasks
from py4j.java_gateway import JavaGateway

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
structure_cpm_dir = "/dat1/nbt2/proj/21-prot/web/data/res/structure_comparison"

################ view ###############
def DUF_SPalign(request):
    jsonFi = "/training/nong/web/data/res/struc_comparison/duf_q2_cov.6.json"
    with open(jsonFi)  as f:
        data_str = json.load(f)
        data = {"data":data_str}
        return JsonResponse(data)
def upload_pdb(request):
    reg_af = re.compile('AF-')
    if request.method == "POST":
        fileDict = request.FILES.items()
        # 获取上传的文件，如果没有文件，则默认为None
        if not fileDict:
            return JsonResponse({'msg': 'No file upload'})

        tempDirName = "file_" + str(uuid.uuid4())
        tempDir = os.path.join(structure_cpm_dir, 'tmp', tempDirName)
        os.system('mkdir -p ' + tempDir)

        for (k, v) in fileDict:
            print("dic[%s]=%s" % (k, v))
            fileData = request.FILES.getlist(k)
            for file in fileData:
                fileName = file._get_name()

                # fileName = reg_W.sub('_', fileName)
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
        
        params = myFunctions.load_POST(request)
        params['dir_1'] = tempDir
        params['proj_type'] = "Structure Comparison"
        res = tasks.structure_comparison.apply_async(args = [params])
        myFunctions.create_submit_form(params, res)
        return JsonResponse({'msg':200})


        # myList = []
        # # alphafold
        # for entry in scanAndFind_notZip(tempDir):
        #     scores = compare_descriptor(entry.path, 'human')
        #     scores2 = compare_descriptor(entry.path, 'pdb')
        #     input_name = entry.name.split('.')[0]
        #     for key, values in scores.items():
        #         lst = defaultdict(str)
        #         lst['input'] = input_name
        #         lst['input_file'] = tempDirName + '/' + entry.name
        #         new_key = key
        #         if reg_af.match(key):
        #             new_key = key.split('-')[1]
        #         lst["target"] = new_key
        #         lst['target_file'] = 'data/alphafold/' + key + ".cif.gz"
        #         lst["geo_score"] = round(values[2], 1)
        #         lst["cn_score"] = round(values[3], 1)
        #         lst["total_score"] = round(values[4], 1)
        #         lst['source'] = 'AlphaFold 2'
        #         myList.append(lst)

        #     # pdb
        #     for key, values in scores2.items():
        #         lst = defaultdict(str)
        #         lst['input'] = input_name
        #         lst['input_file'] = tempDirName + '/' + entry.name
        #         new_key = key
        #         lst["target"] = new_key
        #         lst['target_file'] = 'data/pdb/' + key + ".cif"
        #         lst["geo_score"] = round(values[2], 1)
        #         lst["cn_score"] = round(values[3], 1)
        #         lst["total_score"] = round(values[4], 1)
        #         lst['source'] = 'PDB'
        #         myList.append(lst)

        #     print("match structure AF: ", len(scores))
        #     print("match structure PDB: ", len(scores2))
        # # os.system('rm -rf ')
        # return JsonResponse({'data': myList})

################## function ###############################

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
