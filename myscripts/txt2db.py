#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : 1_format_gene_info.py
# @Date            : 2021/07/26 13:38:41
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, gzip,re



os.environ.setdefault("DJANGO_SETTINGS_MODULE","django_prot.settings")
import django
if django.VERSION >=(1,7):
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
            line = line.strip()
            entry_name, gene_name,  accession , status = line.split('\t')
            info_each = GeneInfo(
                entry_name = entry_name,
                gene_name = gene_name,
                accession = accession,
                status = status
            )
            info.append(info_each)
        f.close()
        GeneInfo.objects.bulk_create(info)

if __name__ == '__main__':
    fire.Fire(RUN)