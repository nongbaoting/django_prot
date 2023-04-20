from django.http import JsonResponse, HttpResponse, FileResponse, Http404
from protein.toolkit import interproscan
from django.core import serializers
import json,os
import pandas as pd
# tada_interproFi = "/dat1/nbt2/proj/23-tadA/work/blast/with_dnaDomain/TadA_all_seqs.fasta.json"
# inter = interproscan.InterproScan()
# tada_interpro = inter.parse(tada_interproFi)
tada_interpro =dict()

def interpro(request):
    global tada_interpro
    if request.method == "GET":
        request_type  = request.GET.get('request_type')
        if request_type == "interpro":
            if  not tada_interpro :
                tada_interproFi = "/dat1/nbt2/proj/23-csr/work/domain/interpro/CSR.fasta_1.json"
                inter = interproscan.InterproScan()
                tada_interpro = inter.parse(tada_interproFi)
            protein_id = request.GET.get('protein_id') 
            print(protein_id)
            return JsonResponse(tada_interpro[protein_id])
        elif request_type == "info_table":
            # table = pd.read_table("/dat1/nbt2/proj/23-tadA/data/results/tadA_interpro-FATCAT_TMalign-withDNA.xls",sep="\t")
            # jsData = df_to_JSjson(table)
            jsData = tableFi_toJSjson("/dat1/nbt2/proj/23-csr/work/domain/interpro/info.tsv")
            # jsData = [item for item in jsData if item["chain_2"]]
            content = get_onePage(request,jsData)
            return JsonResponse(content,safe=False)
        elif request_type == "unidoc":
            f = open("/dat1/nbt2/proj/23-csr/work/domain/unidoc/csr.track.json")
            unidoc = json.loads(json.load(f))
            protein_id = request.GET.get('protein_id') 
            print(protein_id,"unidoc")
            return JsonResponse(unidoc[protein_id],safe=False)
        elif request_type == "Sword2":
            f = open("/dat1/nbt2/proj/23-csr/work/domain/sword/run/null/csr_sword2.track.json")
            unidoc = json.loads(json.load(f))
            protein_id = request.GET.get('protein_id') 
            print(protein_id,"Sword2")
            return JsonResponse(unidoc[protein_id],safe=False)
        else:
            raise Exception("Unknown")
def get_pdbFile(request):
    protein_id = request.GET.get('protein_id') 
    pdbdir= "/dat1/nbt2/proj/23-csr/pdb/colabFoldPdb"
    pdb_fi = os.path.join(pdbdir,f'{protein_id}-cf.pdb')
    if os.path.exists(pdb_fi):
        print(pdb_fi)
        fh = open(pdb_fi, 'rb')
        return FileResponse(fh)
    else:
        print("File not exist!" + pdb_fi)
        raise Http404("File not exist!" + pdb_fi)

import numpy as np
def df_to_JSjson(df):
    columns = df.columns
    print(columns)
    data = []
    for row_index,row in df.iterrows():
        item = {}
        for idx,colname in enumerate(list(columns)):
            value = row[idx]
            if not value or value is np.nan:
                value = 0
            item[colname] = value
        data.append(item)
    return data

def tableFi_toJSjson(fi):
    lines = open(fi).read().strip('\n').split('\n')
    headers = lines[0].split('\t')
    dt = []
    for li in lines[1:]:
        cell = li.split('\t')
        item = {}
        for idx, headname in enumerate(headers):
            value =cell[idx]
            if cell[idx] == 'NA':
                value = ''
            item[headname] = value
        dt.append(item)
    return dt


def get_onePage(request,data):
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))
    field = request.GET.get("field")
    order = request.GET.get('order')
    if field:
        if order == "descending":
            data = sorted(data, key=lambda k: k[field], reverse=True)
        else:
            data = sorted(data, key=lambda k: k[field])
    requestData = data[(currentPage-1) * pageSize: currentPage * pageSize]
    return {"data":requestData, "totalCount": len(data)}