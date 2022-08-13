
from protein.toolkit import TMalign
import os
import pickle
from os import path
from protein.toolkit import *

from django.http import JsonResponse, FileResponse, Http404
struc_cpm_dir = "/dat1/nbt2/proj/21-prot/web/data/res/structure_comparison"
scopeDomain_dir = ''


def getOneItem(request):
    if request.method == "GET":
        myuuid = request.GET.get('uuid')
        dataType = request.GET.get("dataType")
        tmpdir = os.path.join(struc_cpm_dir, myuuid)
        input_pdb = os.path.join(tmpdir, request.GET.get('input_pdb'))
        db_pdb = os.path.join(scopeDomain_dir, request.GET.get('db_pdb'))
        outFi_name = request.GET.get('input_pdb') + request.GET.get('db_pdb')
        item_pickle = os.path.join(tmpdir, outFi_name + '.pickle')
        if dataType == 'info':
            if not os.path.exists(item_pickle):
                item = TMalign.one_to_one(
                    input_pdb, db_pdb, tmpdir, outFi_name)
            else:
                item = myFunctions.pickle_load_file(item_pickle)
            return JsonResponse(item)
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
