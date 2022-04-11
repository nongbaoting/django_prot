import re,uuid,json
import subprocess
from email.mime import base
from os import path
from collections import defaultdict
from protein.models import *

nr_dir = "/dat1/dat/db/nr/divided/all"
base_dir = "/dat1/nbt2/server/PROTsim/django_prot"

def run_phylo(params):
    outdir = create_tmpDir(params['work_dir'], params['uuid'])
    seqFaFi = getNrFasta(params)
    # seqFaFi="/dat1/nbt2/proj/21-prot/web/data/res/phlyo/64376b42-57d7-4d84-a20a-91d48bd27e4e/seqs.fasta"
    sed = "cat seqs.fasta|grep '>'|sed -e 's/>//'|awk '{print $1,$0}' OFS=\"\\t\"|cut -f1 -d']'|awk '{print $0\"]\"}' >id_names.txt;"
    cmd = f"cd {outdir}; mafft --auto --thread 8 {seqFaFi} >seq.msa; \
            iqtree  --fast -T 16 -s seq.msa; \
            {sed} \
            Rscript {base_dir}/myscripts/R/ggtree.R" 
    print(seqFaFi)
    res = subprocess.run(cmd, shell=True)
    status = 0
    return status

def getNrFasta(params):
    tmpDir = create_tmpDir(params['work_dir'], params['uuid'])
    seqIds = params['seqIDs']
    myd = defaultdict(list)
    querySet = NrInfo.objects.filter(protin_id__in=seqIds)
    for q in querySet:
        myd[q.filename].append(q.protin_id)
    
    outFa = path.join(tmpDir, 'seqs.fasta')
    for filename in myd:
        seqListFi = path.join(tmpDir, f'{filename}.fa.list')
        nr_fasta = path.join(nr_dir, filename)
        with open(seqListFi, 'w') as fo:
            fo.write("\n".join(myd[filename]) + '\n')
        cmd = f"seqtk subseq {nr_fasta} {seqListFi} >>{outFa}"
        # print(cmd)
        subprocess.run(cmd, shell=True)
    subprocess.run(f'rm {tmpDir}/*.fa.list', shell=True)
    return outFa

def create_tmpDir(baseDir, myuuid):
    tmpdir = path.join(baseDir, myuuid)
    subprocess.run("mkdir -p " + tmpdir, shell=True)
    return tmpdir