
import prody
from Bio.PDB import *
from protein.toolkit import TMalign, SPalign, Fatcat
import os
import pickle
from os import path
from protein.toolkit import *

from django.http import JsonResponse, FileResponse, Http404
struc_cpm_dir = "/dat1/nbt2/proj/21-prot/web/data/res/structure_alignment/"
database_dir = '/dat1/nbt2/proj/21-prot/dat/pdb/mmCIF/'
struc_result_dir = "/training/nong/protein/res/"

cas9_dir = "/dat1/nbt2/proj/22-cas/work/cas9/dat.1/"
cas9AF_dir = "/training/nong/protein/work/cas9_ColabFold/results/colabFoldPdb/"


def pair_TMalign(request):
    if request.method == "GET":
        dataType = request.GET.get("dataType")
        input_pdb_id = request.GET.get("input_pdb_id")
        target_pdb_id = request.GET.get("target_pdb_id")
        myuuid = input_pdb_id + '_' + target_pdb_id
        print("uuid---", myuuid)
        tmpdir = os.path.join(struc_cpm_dir, myuuid)
        os.system("mkdir -p %s" % tmpdir)
        print(tmpdir)
        input_pdb = os.path.join(cas9_dir, input_pdb_id + '.pdb')
        target_pdb = os.path.join(cas9AF_dir, target_pdb_id + '-cf.pdb')

        outFi_name = f'transform_{myuuid}_tmalign.pdb'
        outFi = os.path.join(tmpdir, outFi_name)
        item_pickle = os.path.join(tmpdir, outFi_name + '.pickle')

        if dataType == 'info':
            if not os.path.exists(item_pickle):
                item = TMalign.pair_align(
                    input_pdb, target_pdb, tmpdir, outFi_name)
            else:
                item = myFunctions.pickle_load_file(item_pickle)
            return JsonResponse(item)
        if dataType == 'input_pdb':
            pdb_fi = os.path.join(tmpdir, outFi_name + '.sup.pdb')
            print(pdb_fi)
            fh = open(pdb_fi, 'rb')
            return FileResponse(fh)
        elif dataType == "db_pdb":
            pdb_fi = target_pdb
            if os.path.exists(pdb_fi):
                print(pdb_fi)
                fh = open(pdb_fi, 'rb')
                return FileResponse(fh)
            else:
                print("File not exist!" + pdb_fi)
                raise Http404("File not exist!" + pdb_fi)


def pair_SPalign(request):
    if request.method == "GET":
        dataType = request.GET.get("dataType")
        input_pdb_id = request.GET.get("input_pdb_id")
        target_pdb_id = request.GET.get("target_pdb_id")
        myuuid = input_pdb_id + '_' + target_pdb_id
        tmpdir = os.path.join(struc_cpm_dir, myuuid)
        os.system("mkdir -p %s" % tmpdir)
        print(tmpdir)
        input_pdb = os.path.join(cas9_dir, input_pdb_id + '.pdb')
        target_pdb = os.path.join(cas9AF_dir, target_pdb_id + '-cf.pdb')

        outFi_name = f'transform_{myuuid}_spalign.pdb'
        outFi = os.path.join(tmpdir, outFi_name)
        item_pickle = os.path.join(tmpdir, outFi_name + '.pickle')

        if dataType == 'info':
            if not os.path.exists(item_pickle):
                spalign = SPalign.pair_align(
                    input_pdb,  target_pdb, tmpdir, outFi_name)
                item = spalign.content

            else:
                item = myFunctions.pickle_load_file(item_pickle)
            return JsonResponse(item)
        if dataType == 'input_pdb':
            pdb_fi = os.path.join(tmpdir, outFi_name)
            print(pdb_fi)
            fh = open(pdb_fi, 'rb')
            return FileResponse(fh)
        elif dataType == "db_pdb":
            pdb_fi = target_pdb
            if os.path.exists(pdb_fi):
                print(pdb_fi)
                fh = open(pdb_fi, 'rb')
                return FileResponse(fh)
            else:
                print("File not exist!" + pdb_fi)
                raise Http404("File not exist!" + pdb_fi)


def pair_Fatcat(request):
    if request.method == "GET":
        dataType = request.GET.get("dataType")
        input_pdb_id = request.GET.get("input_pdb_id")
        target_pdb_id = request.GET.get("target_pdb_id")
        myuuid = input_pdb_id + '_' + target_pdb_id
        tmpdir = os.path.join(struc_cpm_dir, myuuid)
        os.system("mkdir -p %s" % tmpdir)
        print(tmpdir)
        input_pdb = os.path.join(cas9_dir, input_pdb_id + '.pdb')
        target_pdb = os.path.join(cas9AF_dir, target_pdb_id + '-cf.pdb')

        outFi_name = f'transform_{myuuid}_fatcat.pdb'
        outFi = os.path.join(tmpdir, outFi_name)
        item_pickle = os.path.join(tmpdir, outFi_name + '.pickle')
        if dataType == 'info':
            if not os.path.exists(item_pickle):
                fatcat = Fatcat.pair_align(
                    input_pdb,  target_pdb, tmpdir, outFi_name)
                item = fatcat.content

            else:
                item = myFunctions.pickle_load_file(item_pickle)
            return JsonResponse(item)
        if dataType == 'input_pdb':
            pdb_fi = os.path.join(tmpdir, outFi_name)
            print(pdb_fi)
            fh = open(pdb_fi, 'rb')
            return FileResponse(fh)
        elif dataType == "db_pdb":
            pdb_fi = target_pdb
            if os.path.exists(pdb_fi):
                print(pdb_fi)
                fh = open(pdb_fi, 'rb')
                return FileResponse(fh)
            else:
                print("File not exist!" + pdb_fi)
                raise Http404("File not exist!" + pdb_fi)


def select_chain(pdbFi, chain_letter, outdir, db_pdbid):
    outname = os.path.join(outdir, db_pdbid)
    atoms = prody.parseMMCIF(pdbFi, chain=chain_letter)
    outFi = f"{outname}{chain_letter}.pdb"
    if atoms:
        writePrody = prody.writePDB(outFi, atoms)
    return outFi
