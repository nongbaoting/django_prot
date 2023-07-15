
import prody
from Bio.PDB import *
from protein.toolkit import TMalign, SPalign, Fatcat
import os,re
import pickle
from os import path
from protein.toolkit import *

from django.http import JsonResponse, FileResponse, Http404
struc_cpm_dir = "/dat1/nbt2/proj/21-prot/web/data/res/structure_alignment/"
database_dir = '/dat1/nbt2/proj/21-prot/dat/pdb/mmCIF/'
struc_result_dir = "/training/nong/protein/res/"
structure_prediction_dir = "/dat1/nbt2/proj/21-prot/web/data/res/structure_prediction/"
reg_W = re.compile('\W+')

def pair_TMalign(request):
    if request.method == "GET":
        myuuid = request.GET.get('uuid')

        dataType = request.GET.get("dataType")
        tmpdir = os.path.join(struc_cpm_dir, myuuid)
        os.system("mkdir -p %s" % tmpdir)
        print(tmpdir)
        input_pdb = os.path.join(
            struc_result_dir, 'alphafold',  myuuid, 'ranked_0.pdb')
        ## for new structure prediction submit
        if not os.path.exists(input_pdb):
            input_pdb_source = os.path.join(structure_prediction_dir, myuuid, 'AlphaFold2', 'models', 'ranked_0.pdb')
            input_pdb = os.path.join(struc_cpm_dir, 'ranked_0.pdb')
            os.system(f'ln -s {input_pdb_source} {input_pdb}')
            
        db_pdbid = request.GET.get('db_pdbid')
        db_pdbfile = os.path.join(
            database_dir, db_pdbid[1:3].lower(), db_pdbid.lower() + '.cif.gz')

        db_chain = request.GET.get('db_chain')
        db_pdb = select_chain(db_pdbfile, db_chain, tmpdir, db_pdbid)
        outFi_name = f'transform_{myuuid}_{db_pdbid}{db_chain}.pdb'
        outFi = os.path.join(tmpdir, outFi_name)
        item_pickle = os.path.join(tmpdir, outFi_name + '.pickle')
        print(input_pdb, db_pdb)
        if dataType == 'info':
            if not os.path.exists(item_pickle):
                item = TMalign.pair_align(
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


def pair_SPalign(request):
    if request.method == "GET":
        myuuid = request.GET.get('uuid')
        dataType = request.GET.get("dataType")
        myuuid = reg_W.sub('_', myuuid)
        tmpdir = os.path.join(struc_cpm_dir, myuuid)
        os.system("mkdir -p %s" % tmpdir)
        print(tmpdir)
        input_pdb = os.path.join(
            struc_result_dir, 'alphafold',  myuuid, 'ranked_0.pdb')
        ## for new structure prediction submit
        if not os.path.exists(input_pdb):
            input_pdb_source = os.path.join(structure_prediction_dir, myuuid, 'AlphaFold2', 'models', 'ranked_0.pdb')
            input_pdb = os.path.join(struc_cpm_dir, 'ranked_0.pdb')
            os.system(f'ln -s {input_pdb_source} {input_pdb}')

        db_pdbid = request.GET.get('db_pdbid')
        db_pdbfile = os.path.join(
            database_dir, db_pdbid[1:3].lower(), db_pdbid.lower() + '.cif.gz')

        db_chain = request.GET.get('db_chain')
        db_pdb = select_chain(db_pdbfile, db_chain, tmpdir, db_pdbid)
        outFi_name = f'transform_{myuuid}_{db_pdbid}{db_chain}_spalign.pdb'
        item_pickle = os.path.join(tmpdir, outFi_name + '.pickle')
        print(input_pdb, db_pdb)
        if dataType == 'info':
            if not os.path.exists(item_pickle):
                spalign = SPalign.pair_align(
                    input_pdb, db_pdb, tmpdir, outFi_name)
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
            pdb_fi = db_pdb
            if os.path.exists(pdb_fi):
                print(pdb_fi)
                fh = open(pdb_fi, 'rb')
                return FileResponse(fh)
            else:
                print("File not exist!" + pdb_fi)
                raise Http404("File not exist!" + pdb_fi)


def pair_Fatcat(request):
    if request.method == "GET":
        myuuid = request.GET.get('uuid')
        dataType = request.GET.get("dataType")
        tmpdir = os.path.join(struc_cpm_dir, myuuid)
        os.system("mkdir -p %s" % tmpdir)
        print(tmpdir)
        input_pdb = os.path.join(
            struc_result_dir, 'alphafold',  myuuid, 'ranked_0.pdb')

        ## for new structure prediction submit
        if not os.path.exists(input_pdb):
            input_pdb_source = os.path.join(structure_prediction_dir, myuuid, 'AlphaFold2', 'models', 'ranked_0.pdb')
            input_pdb = os.path.join(struc_cpm_dir, 'ranked_0.pdb')
            os.system(f'ln -s {input_pdb_source} {input_pdb}')

        db_pdbid = request.GET.get('db_pdbid')
        db_pdbfile = os.path.join(
            database_dir, db_pdbid[1:3].lower(), db_pdbid.lower() + '.cif.gz')

        db_chain = request.GET.get('db_chain')
        db_pdb = select_chain(db_pdbfile, db_chain, tmpdir, db_pdbid)
        outFi_name = f'transform_{myuuid}_{db_pdbid}{db_chain}_Fatcat.pdb'
        item_pickle = os.path.join(tmpdir, outFi_name + '.pickle')
        print(input_pdb, db_pdb)
        if dataType == 'info':
            if not os.path.exists(item_pickle):
                fatcat = Fatcat.pair_align(
                    input_pdb, db_pdb, tmpdir, outFi_name)
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
            pdb_fi = db_pdb
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
