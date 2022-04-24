

from turtle import end_fill
import django
import os
import sys
import fire
import gzip
import re
from django.core.exceptions import ObjectDoesNotExist
from collections import defaultdict

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "django_prot.settings")
if django.VERSION >= (1, 7):
    django.setup()

class RUN:

    def unique_protCDncbiOne(self,):
        from protein.models import protCDncbiOne
        with open("/dat1/dat/db/nr/CD/ncbi/merge/nr_gt400.1030.merge") as f:
            all = protCDncbiOne.objects.all()
            for li in f:
                protin_id = li.split("\t")[0]
                q = all.filter(protin_id = protin_id)
                if len(q) > 1:
                    q[1].delete()
    
    def unique_protCDncbi(self,):
        from protein.models import protCDncbi
        with open("/dat1/dat/db/nr/CD/ncbi/format/nr_gt400.1030.fa.tsv") as f:
            all = protCDncbi.objects.all()
            for li in f:
                protin_id = li.split("\t")[0]
                q = all.filter(protin_id = protin_id)

                # print(len(q.distinct()))
                myd = defaultdict(str)
                repeat =[]
                for qi in q:
                    # print(qi.id)
                    key = f'{qi.start}-{qi.end}-{qi.cdd_id}'
                    if key in myd:
                        # print(qi)
                        repeat.append(qi)
                    else:
                        myd[key] = qi
                for qrepeat in repeat:
                    qrepeat.delete()
                # break



if __name__ == "__main__":
    fire.Fire(RUN)