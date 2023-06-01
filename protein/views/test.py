from django.http import JsonResponse, HttpResponse, FileResponse, Http404
from protein.toolkit.modules import interproscan
from django.core import serializers
import json,os
def get_pdbFile(request):
    protein_id = request.GET.get('protein_id') 
    pdbdir= "/dat1/nbt2/proj/23-csr/pdb/colabFoldPdb"
    pdb_fi = "/training/nong/protein/work/cas12i/cas12i-4.21/AF-4/alphafold/CN114592040A-cas12Max/CN114592040A-cas12Max/ranked_1.cif"
    if os.path.exists(pdb_fi):
        print(pdb_fi)
        fh = open(pdb_fi, 'rb')
        return FileResponse(fh)
    else:
        print("File not exist!" + pdb_fi)
        raise Http404("File not exist!" + pdb_fi)