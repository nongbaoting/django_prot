#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : json2table.py
# @Date            : 2023/02/03 09:01:09
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip,json
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def blast(seq1, seq2):
    # Create two sequence files
    seq1 = SeqRecord(Seq(seq1),
                    id="seq1")
    seq2 = SeqRecord(Seq(seq2),
                    id="seq2")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")

    # Run BLAST and parse the output as XML
    output = NcbiblastpCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=10)()[0]
    # print(output)
    # blast_result_record = NCBIXML.read(StringIO(output))

    # # Print some information on the result
    # for alignment in blast_result_record.alignments:
    #     for hsp in alignment.hsps:
    #         print(hsp)
    os.remove("seq1.fasta")
    os.remove("seq2.fasta")
    return output.split(",")

class Main:
    def json2table(self, inFi, outFi,inFasta, refFasta):
        faDt = SeqIO.to_dict(SeqIO.parse(inFasta, 'fasta'))
        refrecord = SeqIO.parse(refFasta, 'fasta')
        refsequence = str(next(refrecord).seq)
        # print(refsequence)
        f = open(inFi)
        data = json.loads(json.load(f))

        with open(outFi,'w') as fo:
            fo.write("\t".join(['target', 'target_len', 'desc','identity','align_length']) + '\n')
            for dt in data:
                blast_out = blast(refsequence, str(faDt[dt['target']].seq))
                # print(blast_out)
                fo.write('\t'.join( [dt['target'], str(dt['target_len']),dt['desc'],blast_out[2],blast_out[3] ]) + '\n')
                  
                
if __name__ == '__main__':
    fire.Fire(Main)