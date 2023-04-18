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
#import xmlschema
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
            "trackColor": "#FEFEFE",
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
        track_test ={
            "trackId": 'alternativeTrack2',
            "rowTitle": 'Alt. n°6',
            "trackHeight": 25,
  
            "displayType": 'composite',
            "displayConfig": [
              {
                "displayType": 'block',
                "displayColor": '#27A3B4',
                "displayData": [
                  {
                    "begin": 1,
                    "end": 163,
                    "gaps": [
                      { "begin": 34, "end": 69 },
                      { "begin": 100, "end": 120 },
                    ],
                    "displayId": 'alt_6_domain_0',
                    "sourceId": '#27A3B4',
                    "featureId":
                      'Alt. n°6 • Dom. n°1 | Length: 356 [1-34;212-241;423-714] | Auth. position: 1 - 714',
                  },
                ],
              }
            ]
        }
        test2 = {
            "trackId": 'blockTrack',
            "trackHeight": 20,
            "trackColor": '#FEFEFE',
            "displayType": 'block',
            "displayColor": '#FF0000',
            "rowTitle": 'ECOD Domain',
            "trackData": [
              {
                "featureId": 'a',
                "begin": 1,
                "end": 163,
                "color": '#27A3B4',
               
              },
              
            ],
          },
        rowConfigData.append(seq_track)
        # rowConfigData.append(track_test)
        rowConfigData.append(test2)
        family_track = {
            "trackId": 'Family_track',
            "rowTitle": 'Family',
            "trackHeight": 225,
            "displayType": 'block',
            "trackData": [
            
            ]
        }
        domain_track = {
            "trackId": 'Domain_track',
            "rowTitle": 'DOMAIN',
            
            "displayType": 'block',
            "trackData": [
            
            ]
        }
        homo_superfamily = {
            "trackId": 'HOMOLOGOUS_SUPERFAMILY_track',
            "rowTitle": 'HOMOLOGOUS_SUPERFAMILY',
           
            "displayType": 'block',
            "trackData": [
            
            ]
        }
        conserved_site = []
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
            block_track =  {
                "displayType": 'block',
                "displayColor": '#27A3B4',
                "displayData": [
                    {
                    "begin": start,
                    "end": end,
                    "displayId": sig_acc + str(start) + str(end),
                    "sourceId": '#27A3B4',
                    "featureId":  sig_acc + str(start) + str(end),
                    },
                ],
                }
            block_track =  {
                "featureId": sig_acc,
                "begin": start,
                "end": end,
                "color": '#27A3B4',
                
              },
            if entry:
                entry_acc = entry["accession"]
                entry_name = entry["name"]
                entry_desc = entry["description"]
                entry_type =  entry["type"]
                if entry_type == "FAMILY":
                    family_track["trackData"].append(block_track)
                elif entry_type == "DOMAIN":
                    domain_track["trackData"].append(block_track)
                elif entry_type == "HOMOLOGOUS_SUPERFAMILY":
                    homo_superfamily["trackData"].append(block_track)

        # rowConfigData.extend([family_track,domain_track,homo_superfamily])
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