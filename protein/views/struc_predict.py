from opcode import opname
import re
import os
import time
import json
from protein.models import *
from django.middleware import csrf
from collections import defaultdict
from django.utils import timezone
from django.shortcuts import render
from django.core.files import File
from django.http import JsonResponse, HttpResponse, FileResponse, Http404
from django.core import serializers
import datetime
# Create your views here.

reg_W = re.compile('\W+')
reg_fasta = re.compile('^>\w.*\n[\w]{3,}', re.M)
reg_plainSeq = re.compile('^[A-Za-z]+$', re.M)
reg_genebank = re.compile('^\d+\s+[A-Z]+\s+', re.I)
reg_genebankRep = re.compile('\d+|\s+')
reg_blank = re.compile('\d+|\s+')
uploads_47 = "/training/nong/web/data/uploads/"
uploads_47_struct = "/training/nong/web/data/uploads/structure_predict/"

# struc_result_dir = "/training/nong/protein/db/web/outputs/results/"
struc_result_dir = "/training/nong/protein/res/"


def struc_getFile(request):
    if request.method == "GET":
        filetype = request.GET.get('filetype')
        program = request.GET.get('program')
        job_name = request.GET.get('job_name')
        filename = request.GET.get('filename')
        # TODO新版 job_name 改为uuid
        if filetype == "file":
            filePath = os.path.join(
                struc_result_dir, program,  job_name, filename)
        if filetype == "tar":
            filePath = os.path.join(struc_result_dir, program, filename)
        if os.path.exists(filePath):
            print("file ok: ", filePath)
            fo = open(filePath, 'rb')
            return FileResponse(fo)
        else:
            print("not found: ", filePath)
            raise Http404


def struc_search(request):
    searchType = request.GET.get('searchType')
    searchContent = request.GET.get("searchContent")
    if searchType == 'job_name':
        searchContent = reg_W.sub('_', searchContent)
        queryset = SubmitInfo.objects.filter(job_name__icontains=searchContent)
    elif searchType == 'job_name_id':
        queryset = SubmitInfo.objects.filter(job_name=searchContent)

    elif searchType == 'email_addr':
        queryset = SubmitInfo.objects.filter(email_addr=searchContent)
    elif searchType == 'running_date':
        queryset = SubmitInfo.objects.filter(
            running_date__isnull=False).filter(completed_date__isnull=True)
    elif searchType == 'pending':
        queryset = SubmitInfo.objects.filter(running_date__isnull=True)
    elif searchType == 'completed_date':
        queryset = SubmitInfo.objects.filter(completed_date__isnull=False)

    data = serializers.serialize('json', queryset)
    return JsonResponse(data, safe=False)


re_template = re.compile(
    r'Found an exact template match (?P<templae>\w+).$', re.MULTILINE)


def struc_template(request):

    if request.method == "GET":
        job_name = request.GET.get("job_name")
        program = request.GET.get('program')
        run_log = os.path.join(struc_result_dir, program,
                               job_name, 'run_err.log')
        logs = open(run_log).read()
        groups = re_template.findall(logs)
        data,tmp = [],[]
        for g in groups:
            if g in tmp: continue
            pdbid, chain = g.split('_')
            pdbentry = PDBentry.objects.filter(IDCODE=pdbid).first()
            if pdbentry:
                print(g)
                item = {
                    "CODEID": pdbentry.IDCODE,
                    "Chain":  chain,
                    "HEADER": pdbentry.HEADER,
                    "COMPOUND": pdbentry.COMPOUND
                }
                tmp.append(g)
                data.append(item)
        # data = list(set(data))
        return JsonResponse({"data": data})


def struc_queue(request):
    queryset = SubmitInfo.objects.all().order_by('-upload_date')
    data = serializers.serialize('json', queryset)
    return JsonResponse(data, safe=False)


def structure_submit(request):
    if request.method == 'POST':
        data = {'uploadOk': True}
        body_unicode = request.body.decode('utf-8')
        body = json.loads(body_unicode)
        # TODO: find better way to handle json parsing
        # with open('tmp_alpha.json', 'w') as fo:
        #     fo.write(body_unicode)
        # f = open('tmp_alpha.json')
        # body =json.loads(f.read() )
        # print(body['protein_seq'])
        # print(b2['protein_seq'])
      
        proj_name = body["proj_name"]
        job_name = reg_W.sub('_', proj_name)
        params = {'alphafold': '',
                  'RoseTTAFold_mode': body['RoseTTAFold_mode'],
                  'platform': body['platform']
                  }

        print(job_name)
        # print(request.body.get('protein_seq'))
        seq = body['protein_seq']
        # seq = reg_blank.sub('', seq) + '\n'
        seq = seq.upper()
        print(seq)
        if reg_fasta.match(seq):
            pass
        # elif reg_plainSeq.match(seq):
        #     print('plain')
        #     seq = f'>{job_name}\n' + reg_blank.sub('', seq) + '\n'
        #     print("new_seq:", seq)
        # elif reg_genebank.match(seq):
        #     print('seq is genebank\n', seq)
        #     seq = f'>{job_name}\n' + reg_genebankRep.sub('', seq) + '\n'
        #     print("new_seq:", seq)
        else:
            data = {'uploadOk': False}
            return JsonResponse(data)

        if data['uploadOk']:
            SubmitInfo.objects.create(
                email_addr=body["email"],
                proj_name=body["proj_name"],
                job_name=job_name,
                upload_date=datetime.datetime.now(),
                tools=','.join(body['platform']),
                params=json.dumps(params)
            )
            # writing file into /training/nong/web/public/protein/static/uploads
            upload_fa = os.path.join(uploads_47_struct, f'{job_name}.fa')
            with open(upload_fa, 'w') as f:
                myfile = File(f)
                myfile.write(seq)
                myfile.closed
                print('complete writing files!')

        return JsonResponse(data)
    else:
        context = {'fasta': 'hello world'}
        return JsonResponse(context)


def check_proj(request):
    proj_name = request.GET.get('proj_name')
    pinfo = SubmitInfo.objects.filter(proj_name=proj_name)
    jinof = SubmitInfo.objects.filter(job_name=proj_name)
    data = {"isExist": 0}
    if pinfo or jinof:
        data["isExist"] = 1
    return JsonResponse(data)

# search
def search(request):
    if 'gene' in request.GET:
        # 用不上
        gene = request.GET.get('gene').upper()
        print(gene)
        gene_type = request.GET.get("gene_type")
        q = GeneInfo.objects.filter(gene_name__contains=gene).values()
        # print(q)
        data = []
        for item in q:
            item['index'] = item['gene_name'].index(gene)
            # 应该是有多个的
            a = Structure.objects.filter(accession=item['id'])
            item['alphafold'] = ''
            if a:
                # print(a[0])
                item['alphafold'] = a[0].alphfold
            data.append(item)

        def sort_index(item):
            return item['index']
        data.sort(key=sort_index)

        return JsonResponse({'data': data})
    else:
        return render(request, 'protein/search/search.html', {})


def test(request):
    return render(request, 'protein/test.html', {})


def show_structure(request):
    if 'str_id' in request.GET:
        str_id = request.GET.get("str_id")
        gene_name = request.GET.get('gene_name')
        context = {'str_id': str_id,
                   'gene_name': gene_name,
                   }
    return render(request, 'protein/predict/struc_detail.html', context)


def get_token(request):
    token = csrf.get_token(request)


def index(request):
    context = {'hello': 'hello'}
    return render(request, 'protein/index.html', context)


def alphafold2(request):
    context = {'hello': 'hello'}
    return render(request, 'protein/predict/structure_predict.html', context)


def structure_result(request):
    context = {'hello': 'hello'}
    return render(request, 'protein/predict/struc_detail.html', context)
