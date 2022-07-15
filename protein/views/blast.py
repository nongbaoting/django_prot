
from os import path
import os
import re
import uuid
import time
import json
import ujson

from protein.models import *
from collections import defaultdict
from django.utils import timezone
from django.shortcuts import render
from django.core.files import File
from django.http import JsonResponse, HttpResponse
from django.core import serializers
from django.views.decorators.cache import cache_page
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.core.cache import cache
from protein.toolkit import myFunctions

reg_W = re.compile('\W+')
reg_fasta = re.compile('^>\w.*\n[\w]{3,}', re.M)
reg_plainSeq = re.compile('^[A-Za-z]+$', re.M)
reg_genebank = re.compile('^\d+\s+[A-Z]+\s+', re.I)
reg_genebankRep = re.compile('\d+|\s+')
reg_blank = re.compile('\d+|\s+')

blast_indir = "/dat1/nbt2/proj/21-prot/web/data/uploads/blast"
blast_outdir = "/dat1/nbt2/proj/21-prot/web/data/res/blast/"
blast_outdir_old = "/dat1/nbt2/proj/21-prot/web/data/res/bak/blast"
# watch blast
# from myscripts  import run_blast
# run_blast.watch_dog()

# views
# submit


def psijackhmmer(request):
    data = {}
    status = 200
    if request.method == "POST":
        body_unicode = request.body.decode('utf-8')
        body = json.loads(body_unicode)
        print(body)
        status = 200
        sequence = turn2fasta(body['sequence'])
        body['sequence'] = sequence
        if(not sequence):
            data['msg'] = "Input sequence error!"
            data['status'] = 201
            return JsonResponse(data)
        data['status'] = status
        file_id = str(uuid.uuid4())
        tempFileName = os.path.join(blast_indir, file_id + '.json')
        fafile = os.path.join(blast_indir, file_id + '.fa')
        SubmitInfoNew.objects.create(
            uuid=file_id,
            proj_type="blast",
            email_addr=body["email"],
            job_name=body["job_name"],
            upload_date=datetime.datetime.now(),
            tools=','.join(body['program']),
            params=''
        )

        with open(fafile, 'w') as ff:
            ffa = File(ff)
            ffa.write(sequence)
            # json.dump(body, f)
            # ff.write(sequence)
            ffa.closed

        with open(tempFileName, 'w') as f:
            fjson = File(f)
            fjson.write(body_unicode)
            fjson.closed
        return JsonResponse(data)


blastDict = defaultdict(dict)
blastLt = []


@cache_page(60 * 15)
def res_blast_jackhmmer(request):
    data = {}
    status = 200
    data['status'] = -1
    if request.method == "POST":
        body_unicode = request.body.decode('utf-8')
        body = json.loads(body_unicode)
        print(body)
        myuuid = body["uuid"]
        program = body["program"]
        pageSize = body["pageSize"]
        currentPage = body["currentPage"]
        order = body["order"]
        field = body["field"]
        print("pageSize:", pageSize, currentPage)

        data['msg'] = 'running'
        resFi, res_id = '', ''
        dataset = body['dataset']
        myblast_outdir = blast_outdir
        if dataset == "old":
            oldNew = {
                'cdd_noteCat': 'cdd_notes',
                'cdd_nameCat': 'cdd_annots'

            }
            myblast_outdir = blast_outdir_old
            if field in oldNew:
                field = oldNew[field]

        if program == 'jackhmmer':
            resFi = path.join(myblast_outdir, myuuid,  'jackhmmer.json')
            res_id = '.'.join([myuuid, 'jackhmmer', dataset])
            print(resFi)
        elif program == 'psiblast':
            resFi = path.join(myblast_outdir, myuuid,  'psiblast.json')
            res_id = '.'.join([myuuid, 'psiblast', dataset])
            print(resFi)
        elif program == 'BLASTP':
            print("BLASTP...")
            resFi = path.join(myblast_outdir, myuuid,  'BLASTP.json')
            res_id = '.'.join([myuuid, 'BLASTP', dataset])
            print(resFi)

        dt = []
        if res_id in blastDict:
            dt = blastDict[res_id]
            print(res_id, "in balstDict!")
        elif myuuid != 'upload' and path.exists(resFi):
            print(resFi, "file ok")
            # resFi="/dat1/nbt2/proj/22-cas5/work/blast/out.json"
            f = open(resFi)
            dt = json.loads(json.load(f))
            blastDict[res_id] = dt
            blastLt.append(res_id)
            f.close()
        else:
            "something wrong!"

        print("All:", len(dt))
        if len(blastLt) > 10:
            headKey = blastLt.pop(0)
            print("delte key:", headKey)
            del blastDict[headKey]
        print("keys number: ", len(blastDict.keys()))

        # TODO filter
        if program == 'jackhmmer':
            dt = filter_jackhmmer(body, dt)
        elif program in ['psiblast', "BLASTP"]:
            dt = filter_jackhmmer(body, dt)
        totalCount = len(dt)
        print("After filter: ", totalCount)
        #  order
        print("order field", field)
        if order == "descending":
            dt = sorted(dt, key=lambda k: k[field], reverse=True)
        else:
            dt = sorted(dt, key=lambda k: k[field])

        requestData = dt[(currentPage - 1) * pageSize: currentPage * pageSize]

        if dataset == "old":
            renameData = []
            for item in requestData:
                item['cdd_noteCat'] = item['cdd_notes']
                item['cdd_nameCat'] = item['cdd_annots']
                renameData.append(item)
            requestData = renameData

        content = []
        for items in requestData:
            protin_id = items['target']
            cdd_locs = []
            protin = NrInfo.objects.filter(protin_id=protin_id).first()
            # TODO 旧的用all,新的用NCBI
            regions = protCDncbi.objects.filter(protin_id=protin_id)
            for qr in regions:
                qcd = CDD.objects.get(cdd_id=qr.cdd_id_id)
                cdd_locs.append([qr.start, qr.end, qcd.cdd_name, qcd.cdd_desc])
            items["cdd_locs"] = cdd_locs
            content.append(items)

        content = myFunctions.addCDcorlor(content)
        data = {"totalCount": totalCount,
                "data": content,
                'status': 200
                }

        return JsonResponse(data)


def queue(request):
    pageSize = int(request.GET.get("pageSize"))
    currentPage = int(request.GET.get("currentPage"))
    field = 'job_name'
    if "field" in request.GET:
        field = request.GET.get("field")
        print("field:", field)
        if request.GET.get("order") == "descending":
            field = '-' + field
    info = SubmitInfoNew.objects.filter(
        proj_type="blast").order_by('-upload_date')
    totalCount = info.count()
    # 分页
    paginator = Paginator(info, pageSize)
    print("currentPage, pageSize", currentPage, pageSize)
    try:
        books = paginator.page(currentPage)
    except PageNotAnInteger:
        books = paginator.page(1)
    except EmptyPage:
        books = paginator.page(paginator.num_pages)
    data = serializers.serialize('json', books)

    # requestData = info[(currentPage-1) * pageSize : currentPage * pageSize]
    # data = serializers.serialize('json', requestData)
    data = {"totalCount": totalCount,
            "data": data}
    return JsonResponse(data, safe=False)


def search(request):
    searchType = request.GET.get('searchType')
    searchContent = request.GET.get("searchContent")
    if searchType == 'job_name':
        searchContent = reg_W.sub('_', searchContent)
        queryset = SubmitInfoNew.objects.filter(
            job_name__contains=searchContent)
    elif searchType == 'job_name_id':
        queryset = SubmitInfoNew.objects.filter(job_name=searchContent)
    elif searchType == 'email_addr':
        queryset = SubmitInfoNew.objects.filter(email_addr=searchContent)
    elif searchType == 'running_date':
        queryset = SubmitInfoNew.objects.filter(
            running_date__isnull=False).filter(completed_date__isnull=True)
    elif searchType == 'pending':
        queryset = SubmitInfoNew.objects.filter(running_date__isnull=True)
    elif searchType == 'completed_date':
        queryset = SubmitInfoNew.objects.filter(completed_date__isnull=False)

    data = serializers.serialize('json', queryset)
    data = {"totalCount": queryset.count(),
            "data": data}
    return JsonResponse(data, safe=False)

# 逻辑


def shift_dict(mydt, newKey, lt, keep_len=10):
    lt.append(newKey)
    if len(lt) > keep_len:
        headKey = lt.pop(0)
        del mydt[headKey]
    return [mydt, lt]


def filter_jackhmmer(body, data):
    acc_min = float(body['ident']['min'])
    acc_max = float(body['ident']['max'])
    target_len_min = float(body['target_len']['min'])
    target_len_max = float(body['target_len']['max'])
    print(target_len_max)
    dt = []
    print(body)
    dataset = body['dataset']
    oldNew = {
        'cdd_noteCat': 'cdd_notes',
        'cdd_nameCat': 'cdd_annots'
    }

    for mydt in data:
        if mydt['acc'] >= acc_min and mydt['acc'] <= acc_max and mydt['target_len'] >= target_len_min and mydt['target_len'] <= target_len_max:
            checkDomain = [True]
            for domain in body['domains']:
                if domain['value']:

                    # old
                    if dataset == "old" and domain['type'] == "cdd_noteCat":
                        domain['type'] = 'cdd_notes'
                    if dataset == "old" and domain['type'] == "cdd_nameCat":
                        domain['type'] = 'cdd_annots'

                    isInclude = False
                    dCount = domain['count']
                    keyword = domain['value'].strip()
                    if not domain["exclude"]:
                        # print("key word", keyword, "count: ", dCount,mydt[ domain['type']] )
                        ire_key = f"({keyword}.*?)" + "{" + str(dCount) + ",}"
                        re_k = re.compile(ire_key, re.I)
                        match = re_k.search(mydt[domain['type']])
                        if match:
                            # print('match',match.groups())
                            isInclude = True
                    else:
                        # print('aa')
                        re_exclude = re.compile(keyword, re.I)
                        if not re_exclude.search(mydt[domain['type']]):
                            isInclude = True
                    checkDomain.append(isInclude)
            # print(checkDomain)
            if not False in checkDomain:
                dt.append(mydt)
    return dt


def filter_blast(body, data):
    acc_min = body['ident']['min']
    acc_max = body['ident']['max']
    target_len_min = body['target_len']['min']
    target_len_max = body['target_len']['max']
    dt = []
    for mydt in data:
        # print(type(mydt['acc']))
        if mydt['acc'] > acc_min and mydt['acc'] < acc_max and mydt['target_len'] > target_len_min and mydt['target_len'] < target_len_max:
            for domain in body['domains']:
                if domain['value'] and domain['type'] == 'cdd_name':
                    dCount = domain['count']
                    keyword = domain['value'].strip()
                    if not domain["exclude"]:
                        print("key word", keyword, "count: ", dCount)
                        ire_key = f"({keyword}.*?)" + "{" + str(dCount) + ",}"
                        re_k = re.compile(ire_key, re.I)
                        if re_k.search(mydt[domain['type']]):
                            dt.append(mydt)
                    else:
                        re_exclude = re.compile(keyword, re.I)
                        if not re_exclude.search(mydt[domain['type']]):
                            dt.append(mydt)
    return dt


def checkParams(input):
    if input.find(";") != -1:
        return False
    if input.find('*') != -1:
        return False
    return True


def turn2fasta(seq):
    reg_W = re.compile('\W+')
    reg_fasta = re.compile('^>\w.*\n[\w]{3,}', re.M)
    reg_plainSeq = re.compile('^[A-Za-z]+$', re.M)
    reg_genebank = re.compile('^\d+\s+[A-Z]+\s+', re.I)
    reg_genebankRep = re.compile('\d+|\s+')
    reg_blank = re.compile('\d+|\s+')
    # seq = reg_blank.sub('', seq) + '\n'
    seq = seq.upper()
    print(seq)
    if reg_fasta.match(seq):
        print("seq is fasta")
    elif reg_plainSeq.match(seq):
        print('seq is plain')
        seq = reg_blank.sub('', seq)
        seq = f'>id\n' + reg_blank.sub('', seq) + '\n'

    elif reg_genebank.match(seq):
        print('seq is genebank\n', seq)
        seq = f'>id\n' + reg_genebankRep.sub('', seq) + '\n'
    else:
        return False
    return seq
