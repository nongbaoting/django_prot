from django.test import TestCase

# Create your tests here.

from protein.toolkit import *
import json,fire

class Main:

    def test1(self,paraFile):
        with open(paraFile) as f:
            params = json.loads(json.load(f))
            print(params)
            results = run_pdb_domain_annotations.domain_annotations(params)

if "__name__" == "__main__":
    fire.Fire(Main)