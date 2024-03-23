import re,uuid,json,os
import subprocess
from os import path
from collections import defaultdict
from protein.models import *
from .modules import unidoc,sword2, bindingSites, interproscan, foldseek, binding_CLAPE, metal, idr, PTM
reg_cif = re.compile('cif$')

def domain_annotations_remote( params ):
    '''
    params contains: work_dir,pdbfile,chain
    '''
    work_dir_remote = create_tmpDir(params['pdb_annotations_dir_remote'], params['uuid'])
    target_workdir = os.path.join(params['work_dir'], params['uuid'])

    pdbfile_raw = params['pdbfile']
    pdbfileName = os.path.basename(pdbfile_raw)
  
    target_paramFi = os.path.join(target_workdir, 'params.json' )
    target_pdb = os.path.join(target_workdir,  pdbfileName)
    paramsFi = write_params(params, work_dir_remote)

    # scp paramsFi, pdbfile_raw to 191 server
    os.system(f'scp -p 2042 {paramsFi}    nbt2@s1.z100.vip:{target_paramFi}'  )
    os.system(f'scp -p 2042 {pdbfile_raw} nbt2@s1.z100.vip:{target_pdb}'  )
    ssh_cmd = "ssh -p 2042 nbt2@s1.z100.vip"
    os.system(f"{ssh_cmd} 'cd /dat1/nbt2/server/PROTsim/PISA/django_prot; source ~/.zshrc; conda run -n webProt  python ./test.py test1  {target_paramFi}'")
    # scp json Files to remote server
    os.system(f"scp -p 2042  nbt2@s1.z100.vip:{target_workdir } {work_dir_remote}")

    return 0

def domain_annotations(params):
    '''
    params contains: work_dir,pdbfile,chain
    '''
    work_dir = create_tmpDir(params['work_dir'], params['uuid'])
    pdbfile_raw = params['pdbfile']
    pdbfile = os.path.join( work_dir, 'upload.pdb')
    FaFi = os.path.join( work_dir,'seqs.fasta')
    chain = params['chain']
    write_params(params, work_dir)

    if reg_cif.search(pdbfile_raw):
        pdbfile_change = os.path.join(work_dir,'uploadcif2pdb.pdb')
        cif2pdb(pdbfile_raw, pdbfile_change)
        pdbfile_raw = pdbfile_change
    
    select_chain(pdbfile_raw, chain, pdbfile)
    pdb2fasta(pdbfile, FaFi)

    # run interproscan
    interproscan.run_interproscan(work_dir)
    
    # run unidoc,sword2
    unidoc_track = unidoc.run_unidoc(pdbfile, chain, work_dir)
    sword2_track,sequence = sword2.run_sword2(pdbfile, chain, work_dir)
    
    # run PPI py_scannet
    ppi_track = bindingSites.run_ppi(pdbfile,work_dir)
    dna_track, rna_track, antibody_track = binding_CLAPE.run_clape(pdbfile, work_dir)

    # metal ions
    # change 
    ions = metal.run_LMetalSite(pdbfile,work_dir)

    ## PTM
    ptm = PTM.run_MusiteDeep(pdbfile, work_dir)
     
    ## IDR
    try:
        idr_strack = idr.run_idr(pdbfile,work_dir)
        tracks.append(idr_strack)
    except:
        print("IDR not ok")
        
    # annotation
    foldseek.run_annotate(pdbfile,work_dir)
    
    # tracks
    protvistaTrackFi = os.path.join(work_dir,"protvistData.json")
    tracks = []
    tracks.append(sword2_track['upload'])
    tracks.append(unidoc_track)
    tracks.append(ppi_track)

    tracks.append(dna_track)
    tracks.append(rna_track)
    tracks.append(antibody_track)

    tracks.append(ptm)
    tracks.append(ions)
    
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

def write_params(params,workdir):
    paramsFi = os.path.join(workdir, 'params.json')
    with open(paramsFi, 'w') as fp:
        json.dump(json.dumps(params), fp)
    return paramsFi

def select_chain(inPDB, chain, outPDB):
    cmd = f'pdb_selchain -{chain} {inPDB} > {outPDB}'
    subprocess.run(cmd, shell=True)

def pdb2fasta(pdbfile, FaFi):
    cmd2 = f"pdb2fasta {pdbfile} > {FaFi};sed -i -e '1s/:.*//' {FaFi}"
    subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)
