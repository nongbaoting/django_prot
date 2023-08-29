
# Create your tests here.
import django
import os
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'django_prot.settings')
django.setup()
from protein.toolkit import *
import json,fire

class Main:

    def test1(self,paraFile):
        with open(paraFile) as f:
            params = json.loads(json.load(f))
            print(params)
            results = run_pdb_domain_annotations.domain_annotations(params)


if __name__ == '__main__':
    fire.Fire(Main)