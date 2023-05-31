import re,uuid,json,os
import subprocess
from os import path
from collections import defaultdict
from protein.models import *
from .modules import unidoc,sword2

def domain_annotations(params):
    '''
    params contains: work_dir,pdbfile,chain
    '''
    work_dir = create_tmpDir(params['work_dir'], params['uuid'])
    pdbfile = params['pdbfile']
    chain = params['chain']

    unidoc_track = unidoc.run_unidoc(pdbfile, chain, work_dir)
    sword2_track,sequence = sword2.run_sword2(pdbfile, chain, work_dir)
     
    protvistaTrackFi = os.path.join(work_dir,"protvistData.json")
    tracks = []
    tracks.append(sword2_track['upload'])
    tracks.append(unidoc_track)
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