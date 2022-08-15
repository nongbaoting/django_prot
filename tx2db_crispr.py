#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : tx2db_cripsr.py
# @Date            : 2022/07/15 10:50:27
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os
import sys
import fire
import re
import gzip
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


class Main:

    def CASInfo(self, cas_class, Fi="/dat1/dat/db/uniprot/subset/cas9.txt"):
        from crispr.models import CASInfo
        CASInfo.objects.all().delete()
        info = []
        count, all = 0, 0
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.rstrip('\n').split("\t")
                q = CASInfo(
                    cas_class=cas_class,
                    accession=cell[0],
                    entry_name=cell[1],
                    data_class=cell[2],
                    protein_name=cell[3],
                    gene_name=cell[4],
                    organism=cell[5],
                    taxonomy_id=cell[6],
                    sequence_length=cell[7],
                )
                info.append(q)
                count += 1
                if count > 100000:
                    all += count
                    print("now:", all, "add:", count)
                    CASInfo.objects.bulk_create(info)
                    count = 0
                    info = []
        if count > 0:
            CASInfo.objects.bulk_create(info)
        print(all + count)

    def AlignTMScore(self, Fi="/dat1/nbt2/proj/22-cas/work/cas9/cas9_vs_colafold/result/TMalign/colabFoldPdb.TMalign.txt"):
        '''
        python tx2db_crispr.py  AlignTMScore
        '''
        from crispr.models import AlignTMScore
        AlignTMScore.objects.all().delete()
        info = []
        count, all = 0, 0
        with open(Fi) as f:
            next(f)
            for li in f:
                cell = li.rstrip('\n').split("\t")
                q = AlignTMScore(
                    chain1=cell[0],
                    chain2_acc=cell[1].split('-')[0],
                    chain1_len=int(cell[2]),
                    chain2_len=int(cell[3]),
                    align_len=int(cell[4]),
                    cov_1=float(cell[5]),
                    cov_2=float(cell[6]),
                    RMSD=float(cell[7]),
                    seq_ID=round(float(cell[8]), 2),
                    TMscore_1=round(float(cell[9]), 2),
                    TMscore_2=round(float(cell[10]), 2),
                    d0_1=float(cell[11]),
                    d0_2=float(cell[12]),
                )
                info.append(q)
                count += 1
                if count > 100000:
                    all += count
                    print("now:", all, "add:", count)
                    AlignTMScore.objects.bulk_create(info)
                    count = 0
                    info = []
        if count > 0:
            AlignTMScore.objects.bulk_create(info)
        print(all + count)

    def AlignSPScore(self, tool, delete, Fi="/dat1/nbt2/proj/22-cas/work/cas9/cas9_vs_colafold/result/SPalignNS/colabFoldPdb.SPalignNS_out.txt"):
        '''
        python tx2db_crispr.py AlignSPScore SPalignNS yes
        python tx2db_crispr.py  AlignSPScore SPalign no /dat1/nbt2/proj/22-cas/work/cas9/cas9_vs_colafold/result/SPalign/colabFoldPdb.SPalign_out.txt
        '''
        from crispr.models import AlignSPScore
        if delete == 'yes':
            AlignSPScore.objects.all().delete()
        info = []
        count, all = 0, 0
        with open(Fi) as f:
            headers = next(f).strip().split('\t')
            seqid_idx = headers.index("SEQID")
            for li in f:
                cell = li.rstrip('\n').split("\t")
                q = AlignSPScore(
                    chain1=cell[0],
                    chain2_acc=cell[1].split('-')[0],
                    chain1_len=int(cell[2]),
                    chain2_len=int(cell[3]),
                    SPe=float(cell[4]),
                    SPa=float(cell[5]),
                    SPb=float(cell[6]),
                    tool=tool,
                    eff_len=int(cell[7]),
                    RMSD=float(cell[8]),
                    seq_ID=float(cell[seqid_idx]),
                )
                info.append(q)
                count += 1
                if count > 100000:
                    all += count
                    print("now:", all, "add:", count)
                    AlignSPScore.objects.bulk_create(info)
                    count = 0
                    info = []
        if count > 0:
            AlignSPScore.objects.bulk_create(info)
        print(all + count)

    def AlignFatcatScore(self,
                         Fi="/dat1/nbt2/proj/22-cas/work/cas9/cas9_vs_colafold/result/Fatcat/colabFoldPdb.Fatcat.log"):
        '''
        python tx2db_crispr.py AlignFatcatScore
        '''
        from crispr.models import AlignFatcatScore
        AlignFatcatScore.objects.all().delete()
        info = []
        count, all = 0, 0
        with open(Fi) as f:
            headers = next(f).strip().split('\t')

            for li in f:
                cell = li.rstrip('\n').split("\t")
                q = AlignFatcatScore(
                    chain1=cell[0],
                    chain2_acc=cell[1].split('-')[0],
                    chain1_len=cell[2],
                    chain2_len=cell[3],
                    cov1=cell[4],
                    cov2=cell[5],
                    seq_ID=cell[6],
                    similar=cell[7],
                    alignScore=cell[8],
                    tmScore=cell[9],
                    RMSD=cell[10],
                )
                info.append(q)
                count += 1
                if count > 100000:
                    all += count
                    print("now:", all, "add:", count)
                    AlignFatcatScore.objects.bulk_create(info)
                    count = 0
                    info = []
        if count > 0:
            AlignFatcatScore.objects.bulk_create(info)
        print(all + count)


if __name__ == '__main__':
    fire.Fire(Main)
