#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : 1_format_gene_info.py
# @Date            : 2021/07/26 13:38:41
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

from operator import imod
import django
import os
import sys
import fire
import gzip
import re
from django.core.exceptions import ObjectDoesNotExist

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "django_prot.settings")
if django.VERSION >= (1, 7):
    django.setup()

class RUN:
    def geneinfo(self, Fi="../django_21/sprot_info.txt"):
        from protein.models import GeneInfo
        all = GeneInfo.objects.all()
        all.delete()
        info = []
        f = open(Fi)
        #head = next(f)
        for line in f:
            line = line.strip('\n')
            entry_name, gene_name,  accession, status = line.split('\t')
            info_each = GeneInfo(
                entry_name=entry_name,
                gene_name=gene_name,
                accession=accession,
                status=status
            )
            info.append(info_each)
        f.close()
        GeneInfo.objects.bulk_create(info)

    def Structure(self, Fi="../django_21/sprot_file_relate.txt"):
        from protein.models import Structure  # change
        from protein.models import GeneInfo
        all = Structure.objects.all()  # change
        all.delete()
        info = []
        f = open(Fi)
        # head = next(f)
        for line in f:
            line = line.strip('\n')
            accession, alphfold = line.split('\t')
            q = GeneInfo.objects.get(accession=accession)
            info_each = Structure(  # change
                accession=q,
                alphfold=alphfold
            )
            info.append(info_each)
        f.close()
        Structure.objects.bulk_create(info)  # change

    def AF_uniprot(self, Fi="/dat1/nbt2/proj/21-prot/alphafold/uniprot_info/uniprot-proteome_UP000005640.tab"):
        from protein.models import AF_uniprot
        all = AF_uniprot.objects.all()
        all.delete()

        info = []
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                q = AF_uniprot(
                    uniprot=cell[0],
                    uniprot_name=cell[1],
                    status=cell[2],
                    protein_names=cell[3],
                    gene_names=cell[4],
                    organism=cell[5],
                    length=cell[6]
                )
                info.append(q)

        AF_uniprot.objects.bulk_create(info)

    def ProtCDD(self, Fi="/dat1/dat/db/refseq/refseq.region.cdd.tsv"):
        from protein.models import ProtCDD
        all = ProtCDD.objects.all().delete()
        info = []
        count = 0
        all=0
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                if cell[3]:
                    q =ProtCDD(
                        protin_id=cell[0],
                        length= int(cell[1]),
                        cdd_ids= cell[2],
                        cdd_annots=cell[3],
                        cdd_notes= cell[5],
                        cdd_locs=cell[4],
                    )
                    info.append(q)
                    count +=1
                    if count > 500000:
                        all += count
                        print("now:", all,"add:",count)
                        ProtCDD.objects.bulk_create(info)
                        count=0
                        info=[]
        if count >0:           
            ProtCDD.objects.bulk_create(info)
        print(all + count)
    def ProtCDD2_head(self, Fi="/dat1/dat/db/refseq/head.tsv"):
        from protein.models import ProtCDD2_head
        all = ProtCDD2_head.objects.all().delete()
        info = []
        count = 0
        all=0
        
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                if cell[3]:
                    q =ProtCDD2_head(
                        protin_id=cell[0],
                        length= int(cell[1]),
                        cdd_ids= cell[2],
                        cdd_annots=cell[3],
                        cdd_notes= cell[5],
                        cdd_locs=cell[4],
                    )
                    info.append(q)
                    count +=1
                    if count > 500000:
                        all += count
                        print("now:", all,"add:",count)
                        ProtCDD2_head.objects.bulk_create(info)
                        count=0
                        info=[]
        if count >0:           
            ProtCDD2_head.objects.bulk_create(info)
        print(all + count)
    
    def ProtCDDeach(self, Fi="/dat1/dat/db/refseq/archaea.region_cdd.each.tsv"):
        from protein.models import ProtCDDeach
        all = ProtCDDeach.objects.all().delete()
        info = []
        count = 0
        all=0
        
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                if cell[3]:
                    q =ProtCDDeach(
                        protin_id=cell[0],
                        length= int(cell[1]),
                        cdd_id = cell[2],
                        cdd_name = cell[3],
                        
                        cdd_loca = cell[4],
                        cdd_note = cell[5],
                    )
                    info.append(q)
                    count +=1
                    if count > 500000:
                        all += count
                        print("now:", all,"add:",count)
                        ProtCDDeach.objects.bulk_create(info)
                        count=0
                        info=[]
        if count >0:           
            ProtCDDeach.objects.bulk_create(info)
        print(all + count)
    def CDD(self, Fi = "/dat1/dat/db/CDD/data/cddid.tbl"):
        from protein.models import CDD
        # all = CDD.objects.all().delete()
        info = []
        count = 0
        all=0

        def short(item):
            if len(item ) <50:
                return item
            commers = item.split('. ')
            the_str=commers[0]
            if len(the_str) <10 and len(commers) >=2:
                the_str = '. '.join(commers[0:2])
            if len(the_str) >255:
                the_str = the_str[0:255]
            return the_str
        
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                if cell[3]:
                    q =CDD(
                        pssm =cell[0],
                        cdd_id = cell[1],
                        cdd_name = cell[2],
                        cdd_desc =  cell[3],
                        cdd_desc_short = short(cell[3]),
                        pssm_len = int(cell[4]),
                    )
                    info.append(q)
                    count +=1
                    if count > 500000:
                        all += count
                        print("now:", all,"add:",count)
                        CDD.objects.bulk_create(info)
                        count=0
                        info=[]
        if count >0:           
            CDD.objects.bulk_create(info)
        print(all + count)
    
    def CDDFIx(self,):
        from protein.models import CDD
        Fi = "/dat1/dat/db/CDD/data/cddid.tbl"

        def short(item):
            if len(item ) <50:
                return item
            commers = item.split('. ')
            the_str=commers[0]
            if len(the_str) <10 and len(commers) >=2:
                the_str = '. '.join(commers[0:2])
            if len(the_str) >255:
                the_str = the_str[0:255]
            return the_str

        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                # print(cell[1])
                shor_desc = short(cell[3])
                CDD.objects.filter(cdd_id = cell[1]).update(cdd_desc_short = shor_desc)
      
        
                
    def NrInfo(self, Fi="/dat1/dat/db/nr/divided/nr_gt400_id_files.txt"):
        from protein.models import NrInfo
        nr = NrInfo.objects.all().delete()
        info = []
        count = 0
        all=0
        with open(Fi) as f:
            for li in f:
                cell = li.strip('\n').split("\t")
                q =  NrInfo(
                    protin_id=cell[0],
                    desc= cell[1].split(']')[0] + ']',
                    length= int(cell[2]),
                    filename =  cell[3],
                )
                info.append(q)
                count +=1
                if count > 500000:
                    all += count
                    print("now:", all,"add:",count)
                    NrInfo.objects.bulk_create(info)
                    count=0
                    info=[]
        if count >0:           
            NrInfo.objects.bulk_create(info)
        print(all + count)

    def ProtCDDeach2(self, Fi="/dat1/dat/db/refseq/archaea.region_cdd.each.tsv"):
        from protein.models import ProtCDDeach2
        all = ProtCDDeach2.objects.all().delete()
        info = []
        count = 0
        all=0
        
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                if cell[3]:
                    q =ProtCDDeach2(
                        protin_id=cell[0],
                        length= int(cell[1]),
                        cdd_id = cell[2],
                        cdd_name = cell[3],
                        
                        cdd_loca = cell[4],
                        cdd_note = cell[5],
                    )
                    info.append(q)
                    count +=1
                    if count > 500000:
                        all += count
                        print("now:", all,"add:",count)
                        ProtCDDeach2.objects.bulk_create(info)
                        count=0
                        info=[]
        if count >0:           
            ProtCDDeach2.objects.bulk_create(info)
        print(all + count)
    
    def NrCD(self, Fi="/dat1/dat/db/nr/CD/tmp/blast.tsv"):
        from protein.models import NrCD
        all = NrCD.objects.all().delete()
        info = []
        count = 0
        all=0
        # [query, length, desc, pssm, start, end, evalue,
        #                           biscore, cdd_id, cdd_name, cdd_name_cat, cdd_note_cat
        from protein.models import CDD
        cdd = CDD.objects.all()
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                if cell[3]:
                    cdd_id = cell[8]
                    # print(cdd_id)
                   
                    q =NrCD(
                        cdd_id = cdd.get(cdd_id = cdd_id),
                        protin_id =  cell[0],
                        length =cell[1],
                        desc  =cell[2],
                        pssm  = cell[3],
                        start = int(cell[4]),
                        end = int(cell[5]), 
                        evalue = cell[6],
                        biscore = cell[7],
                        cdd_name = cell[9],
                        cdd_nameCat =cell[10],
                        cdd_idCat = cell[11]
                     
                    )
                    info.append(q)
                    count +=1
                    if count > 100000:
                        all += count
                        print("now:", all,"add:",count)
                        NrCD.objects.bulk_create(info)
                        count=0
                        info=[]
        if count >0:           
            NrCD.objects.bulk_create(info)
        print(all + count)

    def protCDncbi(self, Fi="/dat1/dat/db/nr/CD/tmp/blast.tsv"):
        from protein.models import protCDncbi, NrInfo
        # all = protCDncbi.objects.all().delete() 
        info = []
        count = 0
        all=0
        # [query, length, desc, pssm, start, end, evalue,
        #                           biscore, cdd_id, cdd_name, cdd_name_cat, cdd_note_cat
        from protein.models import CDD
        cdd = CDD.objects.all()
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                if cell[3]:
                    cdd_id = cell[8]
                    # print(cdd_id)
                    q = protCDncbi(
                        cdd_id = cdd.get(cdd_id = cdd_id),
                        cdd_name = cell[9],
                        protin_id = cell[0],
                        length =cell[1],
                      
                        pssm  = cell[3],
                        start = int(cell[4]),
                        end = int(cell[5]), 
                        evalue = cell[6],
                        biscore = cell[7],
                    )
                    info.append(q)
                    count +=1
                    if count > 500000:
                        all += count
                        print("now:", all,"add:",count)
                        protCDncbi.objects.bulk_create(info)
                        count=0
                        info=[]
        if count >0:           
            protCDncbi.objects.bulk_create(info)
        print(all + count)
    
    def protCDncbiCat(self, Fi="/dat1/dat/db/nr/CD/nr_gt400.1030.fa.tsv"):
        from protein.models import protCDncbiCat, NrInfo
        # all = protCDncbiCat.objects.all().delete() 
        info = []
        count,all= 0,0
        nrinfo = protCDncbiCat.objects.all()
        # [query, length, desc, pssm, start, end, evalue,
        #                           biscore, cdd_id, cdd_name, cdd_name_cat, cdd_note_cat
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                q = protCDncbiCat(
                        protin_id =  nrinfo.filter(protin_id = cell[0]).first(),
                        cdd_nameCat =cell[10],
                        cdd_idCat = cell[11]
                )
                info.append(q)
                count +=1
                if count > 100000:
                    all += count
                    print("now:", all,"add:",count)
                    protCDncbiCat.objects.bulk_create(info)
                    count=0
                    info=[]
        if count >0:           
            protCDncbiCat.objects.bulk_create(info)
        print(all + count)

    def protCDncbiOne(self, Fi="/dat1/dat/db/nr/CD/merge/nr_gt400.1030.merge"):
        from protein.models import protCDncbiOne, NrInfo
        # all = protCDncbiOne.objects.all().delete() 
        info = []
        count,all= 0,0
        # [query, length, desc, pssm, start, end, evalue,
        #                           biscore, cdd_id, cdd_name, cdd_name_cat, cdd_note_cat
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.strip('\n').split("\t")
                q = protCDncbiOne(
                        protin_id =  cell[0],
                        protin_nr =  NrInfo.objects.filter(protin_id = cell[0]).first(),
                        cdd_nameCat =cell[1],
                        cdd_idCat = cell[2]
                )
                info.append(q)
                count +=1
                if count > 100000:
                    all += count
                    print("now:", all,"add:",count)
                    protCDncbiOne.objects.bulk_create(info)
                    count=0
                    info=[]
        if count >0:           
            protCDncbiOne.objects.bulk_create(info)
        print(all + count)

    def tempFix(self,):
        from protein.models import CDD,NrCD,protCDncbiOne
        import os
        from os import path
        def scanAndFind_pattern(mydir, mypattern):
            wantFiles = []
            for entry in os.scandir(mydir):
                if entry.is_file()  and mypattern.search(entry.name):
                    wantFiles.append(entry)
                elif entry.is_dir():
                    wantFiles.extend(scanAndFind_pattern(entry.path, mypattern))
            return wantFiles
        import re, json
        re_json = re.compile(r".json$")
        jFiles_dir = "/dat1/nbt2/proj/21-prot/web/data/res/blast/"
        out_dir = "/dat1/nbt2/proj/21-prot/web/data/res/tmp2/"
        jFiles = scanAndFind_pattern(jFiles_dir,re_json)
        fo_l = open("prot_id.txt", 'w')
        lst = open("/dat1/nbt2/proj/21-prot/web/data/res/list.txt", 'r').read().strip('\n').split('\n')
        # from protein.models import ProtCDD
        for entry in jFiles:
            print(entry.path)
            if path.basename(path.dirname(entry.path)) in lst: 
                print("has been Process! continue.....")
                continue
            oFi = path.join(out_dir,  path.basename(path.dirname(entry.path)),entry.name)
            print(oFi)
            os.system("mkdir -p " + path.dirname(oFi))
            with open(entry.path) as f, open(oFi, 'w') as fo:
                myuuid = os.path.basename(os.path.dirname(entry.path))
                # ids = ['6de59a12-69fd-4801-a4e8-a806efd29591',"b86e4ff7-d89c-416d-9ffb-f8a22245e4e7",
                #        'a74a029e-6089-437c-a0ab-999d0a44e547','52153de7-c1a5-4b50-9e2d-00a25cb43a61','3c552ac0-2942-471b-bf2c-43268ef38387']
                # if myuuid in ids:
                dt = json.loads(json.load(f))
                data=[]
                for myd in dt:
                    if int(myd["target_len"] ) >=400:
                        fo_l.write(myd['target'] + "\n")
                    obj = protCDncbiOne.objects.filter(protin_id=myd['target']).first()
                    cdd_names,cdd_ids, cdd_notes = [],[],[]
                    cdd_nameCat, cdd_idCat, cdd_noteCat ='','',''
                    if obj:
                        # print(obj.cdd_annots)
                        cdd_nameCat, cdd_idCat= obj.cdd_nameCat, obj.cdd_idCat
                        # cdd_names = cdd_nameCat.split(', ')
                        cdd_ids = cdd_idCat.split(', ')
                        # print(cdd_nameCat)
                        for cdd_id_ in cdd_ids:
                            q = CDD.objects.get(cdd_id =cdd_id_ )
                            desc = q.cdd_desc_short
                            cdd_notes.append(desc)
                        cdd_noteCat = '^^'.join(cdd_notes)
   
                    myd["cdd_nameCat"] =  cdd_nameCat
                    myd["cdd_idCat"] = cdd_idCat
                    myd["cdd_noteCat"] = cdd_noteCat
                    data.append(myd)   
                dt_json = json.dumps(data)
                json.dump(dt_json, fo)
        fo_l.close()   
if __name__ == '__main__':
    fire.Fire(RUN)
