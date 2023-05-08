# Create your tasks here
from celery import shared_task
import subprocess
from protein.toolkit import *

@shared_task(bind=True, name="AlphaFold2")
def AlphaFold2(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    self.update_state(state='PROGRESS')
    # TODO add funciton here
    
    print("GPU>>>>>>>>")
    # return to results
    return 0

@shared_task(bind=True, name="PDB Domain Annotations")
def pdb_domain_annotations(self, params):
    myuuid = self.request.id
    params['uuid'] = myuuid
    self.update_state(state='PROGRESS')
    # TODO add funciton here
    print(params)
    
    # return to results
    return 0