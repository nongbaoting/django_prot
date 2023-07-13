from .modules import structure_prediction
import os
from  protein.toolkit import myFunctions

def predict(params):
    '''
    params contains: work_dir,pdbfile,chain
    '''
    work_dir = myFunctions.create_tmpDir(params['work_dir'], params['uuid'])
    protein_fasta = os.path.join(work_dir, 'protein.fasta')
    nucleic_fasta = os.path.join(work_dir, 'nucleic.fasta')
    nucleic_type= params['nucleic_type']
    project_type = params['project_type']
    with open(protein_fasta, 'w') as f_protein, open(nucleic_fasta, 'w') as f_nucleic:
        f_protein.write(params['protein_seq'])
        f_nucleic.write(params['nucleic_seq'])
    if project_type == "RNP":
        rf2na = structure_prediction.RoseTTAFold2NA(work_dir)
        rf2na.run(protein_fasta, nucleic_type, nucleic_fasta)
    return 0