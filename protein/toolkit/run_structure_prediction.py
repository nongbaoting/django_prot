from .modules import structure_prediction
import os,json
from  protein.toolkit import myFunctions

def predict(params):
    '''
    params contains: work_dir,pdbfile,chain
    '''
    work_dir = myFunctions.create_tmpDir(params['work_dir'], params['uuid'])
    with open(f'{work_dir}/params.json', 'w') as fo:
        json.dump(json.dumps(params), fo)

    protein_fasta = os.path.join(work_dir, 'protein.fasta')
    project_type = params['project_type']
    with open(protein_fasta, 'w') as f_protein:
        f_protein.write(params['protein_seq'])
        
    if project_type == "RNP":
        nucleic_fasta = os.path.join(work_dir, 'nucleic.fasta')
        with open(nucleic_fasta, 'w') as f_nucleic:
            nucleic_type= params['nucleic_type']
            f_nucleic.write(params['nucleic_seq'])
        rf2na = structure_prediction.RoseTTAFold2NA(work_dir, params)
        rf2na.run(protein_fasta, nucleic_type, nucleic_fasta)
    elif project_type == "Protein":
        alphafold2 = structure_prediction.AlphaFold2(work_dir, params)
        alphafold2.run(protein_fasta)
    else:
        return 1
    return 0