import json
from collections import defaultdict

def alignFatcatScore(request):
    # 分页
    pageSize = int(request.GET.get('pageSize'))
    currentPage = int(request.GET.get('currentPage'))

    tool = request.GET.get('tool')
    field = '-seq_ID'
    filters = json.loads(request.GET.get('filters'))
    print(filters)
    if request.GET.get("field"):
        print(request.GET.get('field'))
        field = re_field.sub('', request.GET.get("field"))
        print("field:", field)
        if request.GET.get("order") == "descending":
            field = '-' + field

    objs = AlignFatcatScore.objects.all().filter(
        chain1__in=cas9Dt[filters['protein']],).all().filter(chain2_len__gt=filters['min_len'],
                                                         seq_ID__gt=filters['min_SI'], seq_ID__lt=filters['max_SI'],
                                                         chain2_len__lt=filters['max_len']).order_by(field)

    if filters['exclude_knownCas'] == True:
        objs = filter_known_cas(objs)

    if filters['candidates'] == True:
        objs = objs.filter(chain2_acc__in=candidates)
    totalCount = objs.count()
    print('-----------', totalCount)
    requestData = objs[(currentPage-1) * pageSize: currentPage * pageSize]
    data = []
    for item in requestData:
        it = model_to_dict(item)
        # print(it)
        cas9 = CASInfo.objects.get(accession=item.chain2_acc)
        it["protein_name"] = cas9.protein_name
        it["organism"] = cas9.organism
        it["genome_genbank"] = cas9.genome_genbank
        it["protein_genebankID"] = cas9.protein_genebankID
        it['repeatinfo'] = repeatDt[cas9.genome_genbank]
        data.append(it)
    # data = serializers.serialize('json', requestData)
    content = {"totalCount": totalCount,
               "data": data}
    return JsonResponse(content)

candidates = []
def parse_info(infoFile = "/dat1/nbt2/proj/22-cas/work/cas12/comparison/result/Fatcat/colabFoldPdb.Fatcat.PAM.xls"):
    data = []
    with open(infoFile, 'r') as f:
        header = next(f).strip('\n').split("\t")
        header_len = len(header)
        for li in f:
            cell = li.split('\n').split("\t")
            chain1 = cell[0]
            chain2_acc = cell[0]
            chain1_len = cell[0]
            chain2_len = cell[0]
            cov1 =cell[0]
            cov2 =cell[0]
            seq_ID =cell[0]
            similar =cell[0]
            alignScore =cell[0]
            tmScore =cell[0]
            RMSD =cell[0]
            dt = {

            }

        data.append(dt)
    return data
