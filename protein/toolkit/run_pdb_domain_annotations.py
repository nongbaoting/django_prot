import re,uuid,json
import subprocess
from os import path
from collections import defaultdict
from protein.models import *
from .modules import unidoc,sword2

def domain_annotations(params):
    work_dir = create_tmpDir(params['work_dir'], params['uuid'])
    pdbfile = params['pdbfile']
    chain = params['chain']

    status = unidoc.run_unidoc(pdbfile, chain, work_dir)
    sword2.run_sword2(pdbfile, chain, work_dir)

    return status





####
def create_tmpDir(baseDir, myuuid):
    tmpdir = path.join(baseDir, myuuid)
    subprocess.run("mkdir -p " + tmpdir, shell=True)
    return tmpdir