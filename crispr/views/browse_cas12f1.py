import json,re
from collections import defaultdict
from django.http import JsonResponse, FileResponse, Http404

candidates = []
infoFile2 = "/dat1/nbt2/proj/22-cas/work/cas12/comparison/result/Fatcat/colabFoldPdb.Fatcat.PAM.xls"
infoFile ="/dat1/nbt2/proj/22-cas/work/cas12/comparison_1027/pam/results/cas12f-1027_tracer_info_pam.xls"
infoFile_asCluster = "/dat1/nbt2/proj/22-cas/work/cas12/run221104_genbank/length_350_450/clusterWithAscas12f/comparison/result/Fatcat/colabFoldPdb.Fatcat.PAM.xls"
def parse_info(infoFile = infoFile):
    data = []
    with open(infoFile, 'r') as f:
        header = next(f).strip('\n').split("\t")
        header_len = len(header)
        n=0
        if 'NA' in header:
            n=1
        for li in f:
            cell = li.strip('\n').split("\t")
            
            dt = {
            "chain1": cell[0],
            "chain2_acc": cell[1].split("-cf")[0],
            "chain1_len": float(cell[2]),
            "chain2_len": float(cell[3]),
            "cov1":float(cell[4]),
            "cov2":float(cell[5]),
            "seq_ID":float(cell[6]),
            "similar":float(cell[7]),
            "alignScore":float(cell[8]),
            "tmScore":float(cell[9]),
            "RMSD":float(cell[10]),
            "taxid":cell[11+n],
            "organism":cell[12+n],
            "upPam":cell[13+n],
            "upScore": round(float(cell[14+n]),2),
            "downPam":cell[15+n],
            "downScore": round(float(cell[16+n]),2)
            }
            if n==1:
                dt['tracer'] = cell[11]
            else:
                dt['tracer'] ='not known'
            data.append(dt)
    return data

cas12f1 = parse_info()
cas12f1_2 = parse_info(infoFile2)

cas12f1_3 = parse_info(infoFile_asCluster)
cas12f1.extend(cas12f1_2)
cas12f1.extend(cas12f1_3)
re_field = re.compile(r'fields.')

def alignFatcatScore(request):
    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))

    tool = request.GET.get('tool')
    field = '-seq_ID'
    filters = json.loads(request.GET.get('filters'))
    print(filters)
    cas12f1_data = cas12f1
    if request.GET.get("field"):
        print(request.GET.get('field'))
        field = re_field.sub('', request.GET.get("field"))
        print("field:", field)
        if request.GET.get("order") == "descending":
            cas12f1_data = sorted(cas12f1, key=lambda x:x[field], reverse=True)
        else:
            cas12f1_data = sorted(cas12f1, key=lambda x:x[field])
    cas12f1_filter = []
    for item in cas12f1_data:
        if item['chain1' ] == filters['protein']:
            if item['chain2_len'] >= float(filters['min_len']) and item['chain2_len'] <= float(filters['max_len']):
                if item['seq_ID'] >= float(filters['min_SI'] )and item['seq_ID']  <= float(filters['max_SI']):
                    if item['organism']:
                        if re.search(filters['organism'],  item['organism']):
                            cas12f1_filter.append(item)
                    else:
                        cas12f1_filter.append(item)
    # objs = AlignFatcatScore.objects.all().filter(
    #     chain1__in=cas9Dt[filters['protein']],).all().filter(chain2_len__gt=filters['min_len'],
    #                                                      seq_ID__gt=filters['min_SI'], seq_ID__lt=filters['max_SI'],
    #                                                      chain2_len__lt=filters['max_len']).order_by(field)

    # if filters['exclude_knownCas'] == True:
    #     objs = filter_known_cas(objs)

    # if filters['candidates'] == True:
    #     objs = objs.filter(chain2_acc__in=candidates)
    totalCount = len(cas12f1_filter)
    print('-----------', totalCount)
    requestData = cas12f1_filter[(currentPage-1) * pageSize: currentPage * pageSize]
    # data = []
    # for item in requestData:
    #     it = model_to_dict(item)
    #     # print(it)
    #     cas9 = CASInfo.objects.get(accession=item.chain2_acc)
    #     it["protein_name"] = cas9.protein_name
    #     it["organism"] = cas9.organism
    #     it["genome_genbank"] = cas9.genome_genbank
    #     it["protein_genebankID"] = cas9.protein_genebankID
    #     it['repeatinfo'] = repeatDt[cas9.genome_genbank]
    #     data.append(it)
    # data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": requestData}
    return JsonResponse(content)

