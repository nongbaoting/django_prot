#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : PDB.py
# @Date            : 2022/01/08 16:43:05
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os,subprocess
import sys
import fire
import re
from Bio.PDB import PDBParser, MMCIFIO
from collections import defaultdict
from pdbecif.mmcif_tools import MMCIF2Dict
from pdbecif.mmcif_io import CifFileReader, CifFileWriter
# import pdbecif.mmcif_io as mmcif

reg_blank = re.compile(r'\s+')
pattern_pdb = re.compile('.pdb$|.pdb.gz$', re.IGNORECASE)
pattern_cif = re.compile('.cif$|cif.gz$', re.IGNORECASE)

def scanAndFind_pattern(mydir, mypattern):
    wantFiles = []
    for entry in os.scandir(mydir):
        if (entry.is_file() ) and mypattern.search(entry.name):
            wantFiles.append(entry)
        elif entry.is_dir():
            wantFiles.extend(scanAndFind_pattern(entry.path, mypattern))
    return wantFiles

class Main:

    def AF_pdb2cif(self,input_pdb,output_cif):
        
        parser = PDBParser()

        # Parse the PDB file
        structure = parser.get_structure("my_structure", input_pdb)

        # Create an mmCIF output object
        mmcif_writer = MMCIFIO()

        # Save the structure in mmCIF format
        mmcif_file = input_pdb + '.cif'
        mmcif_writer.set_structure(structure)
        mmcif_writer.save(mmcif_file)


        
        mmcif_dict = MMCIF2Dict()
        cif_dict = mmcif_dict.parse(mmcif_file)

        
        chains = sorted(list(set(cif_dict['my_structure']['_atom_site']['label_asym_id'])))
        chains_entity = {i:idx +1 for idx,i in enumerate(chains)}

        label_asym_id=[]
        label_comp_id=[]
        label_seq_id=[]
        metric_id=[]
        metric_value=[]
        model_id=[]
        ordinal_id=[]
        label_entity_id = 0
        label_asym_id_each = ''
        for idx, asym_id in enumerate(cif_dict['my_structure']['_atom_site']['label_asym_id']):
            entity_id = chains_entity[asym_id]
            if cif_dict['my_structure']['_atom_site']['label_asym_id'][idx] != label_asym_id_each:
                label_asym_id_each =  cif_dict['my_structure']['_atom_site']['label_asym_id'][idx]
                label_entity_id +=1
            # cif_dict['my_structure']['_atom_site']['label_entity_id'][idx] = str(label_entity_id)
            cif_dict['my_structure']['_atom_site']['label_entity_id'][idx] = 1
            label_asym_id.append(cif_dict['my_structure']['_atom_site']['label_asym_id'][idx])
            label_comp_id.append(cif_dict['my_structure']['_atom_site']['label_comp_id'][idx])
            label_seq_id.append(cif_dict['my_structure']['_atom_site']['label_seq_id'][idx])
            metric_id.append(2)
            metric_value.append( '%.2f' % ( float( cif_dict['my_structure']['_atom_site']['B_iso_or_equiv'][idx])))
            model_id.append(1)
            ordinal_id.append(cif_dict['my_structure']['_atom_site']['label_seq_id'][idx])
        # max_res = max([int(i) for i in label_seq_id])
        # res_occupy = len(str(max_res))
        # label_seq_id_left = [ "{:<{}}".format(i, res_occupy) for i in label_seq_id]
        cif_dict['my_structure']['_ma_qa_metric_local'] = {
                "label_asym_id":label_asym_id,
                "label_comp_id":label_comp_id,
                "label_seq_id":label_seq_id,
                "metric_id":metric_id,
                "metric_value":metric_value,
                "model_id":model_id,
                "ordinal_id":label_seq_id,
            }
        cif_dict['my_structure']['_ma_qa_metric'] ={
            'id':[1,2],
            'mode':['global','local'],
            'name':['pLDDT','pLDDT'],
            'software_group_id':[1,1],
            'type':['pLDDT','pLDDT'],
        }
        cfd1 = CifFileWriter(output_cif)
        cfd1.write(cif_dict)
        os.remove(mmcif_file)

    def scanAF_pdb2cif(self,mydir):
        for entry in scanAndFind_pattern(mydir,pattern_pdb):
            output_cif = os.path.join( os.path.dirname(entry.path), entry.name.split('.pdb')[0]) + '.cif'
            print(entry.path,'-->',output_cif)
            self.AF_pdb2cif(entry.path, output_cif)

    def scanAndSelectRegion(self,region,mydir,outdir,chain='A'):
        start,end = region
        for entry in scanAndFind_pattern(mydir,pattern_pdb):
            output_pdb = os.path.join( outdir, entry.name.split('.pdb')[0]) + f"-{start}_{end}"+ '.pdb'
            cmd = f"pdb_selchain -A {entry.path}|pdb_selres -{start}:{end} > {output_pdb}"
            print(cmd)
            subprocess.run(cmd, shell=True)
    def membrancePDBAddChain(self,inputPDB, outputPDB):
        pass
          

if __name__ == '__main__':
    fire.Fire(Main)
