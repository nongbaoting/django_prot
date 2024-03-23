# Create your tasks here
from celery import shared_task
import subprocess
from protein.toolkit import *
from protein.toolkit.modules import PISA_TMalign
@shared_task(bind=True, name="Phylogenetic Analysis",property=9)
def run_phylo(self, params):
    myuuid = self.request.id
    self.update_state(state='PROGRESS')
    params['uuid'] = myuuid
    status = run_phylogenetic.run_phylo(params)
    return status

@shared_task(bind=True, name="Download Fasta")
def getNrFasta(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    status = run_phylogenetic.getNrFasta(params)
    return status
    
@shared_task(bind=True, name="CD Search",property=0)
def cd_search_save(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    self.update_state(state='PROGRESS')
    status = run_cdSearch.search_save(params)
    # return to results
    return status

@shared_task(bind=True, name="Structure Comparison")
def structure_comparison(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    self.update_state(state='PROGRESS')
   
    tmalgin = TMalign.TMalgin(params)
    result = tmalgin.run()
    print(result)
    # return to results
    return result

@shared_task(bind=True, name="PDB Domain Annotations")
def pdb_domain_annotations(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    self.update_state(state='PROGRESS')
   
    print(params)
    results =''
    results = run_pdb_domain_annotations.domain_annotations(params)
    # if params['ip'] == "172.22.148.191":
    #     results = run_pdb_domain_annotations.domain_annotations(params)
    # else:
    #     results = run_pdb_domain_annotations.domain_annotations_remote(params)
    return results

from myscripts import run_blast 

@shared_task(bind=True, name="blast")
def blast(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    self.update_state(state='PROGRESS')
    
    print(params)
    results = run_blast.run_blast(params)
    return results

@shared_task(bind=True, name="Structure prediction", concurrency=1)
def structure_prediction(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    self.update_state(state='PROGRESS')
    print("uuid: ", myuuid)
    print(params)
    results = run_structure_prediction.predict(params)
    return results

@shared_task(bind=True, name="TMalign", concurrency=1)
def Alignment(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    self.update_state(state='PROGRESS')
    print("uuid: ", myuuid)
    print(params)
    results = PISA_TMalign.run_TMalign(params)
    return results