#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : format.py
# @Date            : 2023/04/17 14:46:38
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip,json,pickle
from collections import defaultdict
# from protein.toolkit.myFunctions import colorSet3D
#import xmlschema

colorDt = {
"CDD":'#8DD3C7',
"FUNFAM":'#FCCDE5',
"GENE3D":"#27A3B4",
"PANTHER":"#1B9E77",
"PFAM":'#FCCDE5',
"PHOBIUS":'#FFFFB3',
"PIRSR":'#BEBADA',
"PROSITE_PATTERNS":'#D9D9D9',
"PROSITE_PROFILES":'#D9D9D9',
"SFLD" :"#1B9E77",
"SMART":'#D95F02',
"SUPERFAMILY":'#66A61E',
"TIGRFAM":'#666666',
"TMHMM":'#8C510A',
}

class InterproScan:
    def __init__(self):
        self.dt = dict()
       
    def parse(self, jsonFi):
        with open(jsonFi, "r") as f:
            self.dt = json.loads(f.read())
            # print(self.dt)
        return self.interpro2rcsb()
    def interpro2rcsb(self):
        results =self.dt["results"]
        item = results[0]
        sequence = item["sequence"]
        matches = item["matches"]
        rowConfigData =[]
        seq_track =  {
            "trackId": "sequenceTrack",
            "trackHeight": 20,
            # "trackColor": "#FEFEFE",
            "displayType": "sequence",
            "nonEmptyDisplay": True,
            "rowTitle": "SEQUENCE",
            "trackData": [
                {
                    "begin": 1,
                    "value": sequence,
                },
                ],
            }
       
        rowConfigData.append(seq_track)
        family_track = {
            "trackId": "Family_track",
            "rowTitle": "Family",
            "trackHeight": 25,
            "titleFlagColor":"#D83E7C",
            # "borderWidth": 10,
            "displayColor":"#D83E7C",
            "displayType": "block",
            "trackData": [
            
            ]
        }
        domain_track = {
            "trackId": "Domain_track",
            "rowTitle": "DOMAIN",
            "displayType": "block",
            "titleFlagColor":"#27A3B4",
            "trackData": [
            
            ]
        }
        homo_superfamily = {
            "trackId": "HOMOLOGOUS_SUPERFAMILY_track",
            "rowTitle": "SUPERFAMILY",
            "displayType": "block",
            "titleFlagColor":'#C08423',
            "trackData": [
            
            ]
        }
        CONSERVED_SITE = {
            "trackId": "CONSERVED_SITE_track",
            "rowTitle": "CONSERVED_SITE",
            "displayType": "block",
            "titleFlagColor":'#C08423',
            "trackData": [
            
            ]
        }
        conserved_site = []
        tracks_num = 0
        for match in matches:
            signature = match["signature"]
            locations = match["locations"]
            model_ac = match["model-ac"]
            # evalue = match["evalue"]
            sig_acc = signature["accession"]
            sig_name = signature["name"]
            sig_desc = signature["description"]
            sig_source = signature["signatureLibraryRelease"]["library"]
            locations = match["locations"]
            start,end = locations[0]["start"],locations[0]["end"]
            entry = signature["entry"]
            tracks_num +=1
            # print(sig_source)
            block_track =  {
                "featureId": f'{sig_source} | {sig_acc} | {sig_desc}',
                "begin": start,
                "end": end,
                "color": colorDt[sig_source],
              }
            if entry:
                entry_acc = entry["accession"]
                entry_name = entry["name"]
                entry_desc = entry["description"]
                entry_type =  entry["type"]
                print(entry_type)
                if entry_type == "FAMILY":
                    family_track["trackData"].append(block_track)
                elif entry_type == "DOMAIN":
                    domain_track["trackData"].append(block_track)
                elif entry_type == "HOMOLOGOUS_SUPERFAMILY":
                    homo_superfamily["trackData"].append(block_track)
                elif entry_type == "CONSERVED_SITE":
                    CONSERVED_SITE["trackData"].append(block_track)


        rowConfigData.extend([homo_superfamily,family_track,domain_track,CONSERVED_SITE])
        data= {
            "sequence": sequence,
            "rowConfigData": rowConfigData
            }
        return data

class Main:
    def interproscan(self,fi="CDKAL_HUMAN.json"):
        inter = InterproScan()
        rowConfigData = inter.parse(fi)
       
       

if __name__ == "__main__":
    fire.Fire(Main)