import re,uuid,json,os
import subprocess
from os import path
from collections import defaultdict
from protein.models import *
from .modules import unidoc,sword2,bindingSites,interproscan
reg_cif = re.compile('cif$')
def domain_annotations(params):
    '''
    params contains: work_dir,pdbfile,chain
    '''
    work_dir = create_tmpDir(params['work_dir'], params['uuid'])
    pdbfile = params['pdbfile']
    chain = params['chain']

    if reg_cif.search(pdbfile):
        pdbfile_change = os.path.join(work_dir,'upload.pdb')
        cif2pdb(pdbfile, pdbfile_change)
        pdbfile = pdbfile_change

    # run unidoc,sword2
    unidoc_track = unidoc.run_unidoc(pdbfile, chain, work_dir)
    sword2_track,sequence = sword2.run_sword2(pdbfile, chain, work_dir)
    
    # run PPI py_scannet
    ppi_track = bindingSites.run_ppi(pdbfile,work_dir)

    # run interproscan
    interproscan.run_interproscan(work_dir)
    
    protvistaTrackFi = os.path.join(work_dir,"protvistData.json")
    tracks = []
    tracks.append(sword2_track['upload'])
    tracks.append(unidoc_track)
    tracks.append(ppi_track)
    protvistaData = {
        "displayNavigation": True,
        "displaySequence": True,
        "sequence": sequence,
        "length": len(sequence), 
        "offset": 1,
        "tracks":tracks,
    }
    with open(protvistaTrackFi, 'w') as fp:
        json.dump(json.dumps(protvistaData), fp)
    return 0



####
def create_tmpDir(baseDir, myuuid):
    tmpdir = path.join(baseDir, myuuid)
    subprocess.run("mkdir -p " + tmpdir, shell=True)
    return tmpdir

## cif2pdb
import Bio.PDB
def cif2pdb(cifFile,pdbFile):
    parser = Bio.PDB.MMCIFParser()
    structure = parser.get_structure("name", cifFile)
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdbFile)