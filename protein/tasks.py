# Create your tasks here
from __future__ import absolute_import, unicode_literals
from celery import shared_task
import subprocess
from protein.toolkit import *



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
    # TODO add funciton here
    tmalgin = TMalign.TMalgin()
    # return to results
    return 0
