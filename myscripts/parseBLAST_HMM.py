
from Bio.Blast import NCBIXML
from protein.models import ProtCDD
import os, sys, fire, re,json
from collections import defaultdict
reg_sp = re.compile("\s+")
reg_cm = re.compile("#")
reg_CRISPR = re.compile("CRISPR", re.I)
from protein.models import CDD, protCDOne

def retrieve_cdd(protin_id):
    obj = protCDOne.objects.filter(protin_id=protin_id).first()
    cdd_names,cdd_ids, cdd_notes = [],[],[]
    cdd_nameCat, cdd_idCat, cdd_noteCat ='','',''
    if obj:
        cdd_nameCat, cdd_idCat= obj.cdd_nameCat, obj.cdd_idCat
        cdd_ids = cdd_idCat.split(', ')
        for cdd_id_ in cdd_ids:
            q = CDD.objects.get(cdd_id =cdd_id_ )
            desc = q.cdd_desc_short
            # cdd_notes.append(desc)
            cdd_names.append(q.cdd_name)
        # cdd_noteCat = '^^'.join(cdd_notes)
        cdd_nameCat = ', '.join(cdd_names)
    return [ cdd_nameCat, cdd_idCat, cdd_noteCat]

def parse_hmmer(Fi):
    myd = defaultdict(list)
    blocks = open(Fi, 'r').read().strip().split(">>")
    # print(blocks[1])
    for blk in blocks[1:]:
        seq_id = blk.split("\n")[0].strip().split(" ")[0].strip()
        myd[seq_id].append(blk) 
    return myd

def parse_hmmerdomtbl(Fi):
    tbl = defaultdict(dict)
    sameDomain =defaultdict(list)
    seqWithAnnotate = 0
    with open(Fi, 'r') as f:
        for li in f:
            li = li.strip()
            if re.match("#", li):continue
            cell = reg_sp.split( li)
            target, _, tlen,query, _, qlen,E_value, score,bias, domain,domain_numbers, c_Evalue, i_Evalue, domain_score, domain_bias, query_from, query_to, target_from, target_to, env_from, env_to, acc = cell[0:22]
            tlen = int(tlen)
            # if tlen <= 400 : continue
            new_cell = cell[0:22].copy()
            desc = li[181:]
            cdd_names, cdd_ids, cdd_notes = retrieve_cdd(target )
            dt = {
                "target": target,
                "query": query,
                "target_len": tlen,
                "E_value": float(E_value),
                "score": float(score),
                "domain":domain,
                "domain_count": domain_numbers,
                "acc": float(acc),
                'desc': desc,
                "cdd_nameCat" :cdd_names ,
                "cdd_noteCat" :cdd_notes,
                "cdd_idCat" : cdd_ids,
            }
            
            sameDomain[cdd_names].append([target, tlen, desc])
            tbl[target] = dt
    newSameDomainList = []
    totalSequence = len(tbl.keys())

    for cdd_nameCat in sameDomain:
        min_tlen, max_tlen = 10000000,0
        descDt = defaultdict(int)
        proteinLenDt = defaultdict(list)
        count = len(sameDomain[cdd_nameCat])
        for target, tlen, desc in sameDomain[cdd_nameCat]:
            if cdd_nameCat != "":
                seqWithAnnotate +=1
            proteinLenDt[str(tlen)].append(target)
            min_tlen = min(tlen,min_tlen)
            max_tlen =max(tlen,max_tlen)
            descDt[desc.split('[')[0].split(',')[0] ] +=1
        commonDesc = sorted(descDt.items(), key=lambda k: k[1], reverse=True)[0][0]
        lenMode = sorted(proteinLenDt.keys(), key=lambda k: int(k), reverse=True)[0]
        rep_protein_id = proteinLenDt[lenMode][0]
        domainInfo = {
            "cdd_nameCat": cdd_nameCat,
            "min_tlen": min_tlen,
            "max_tlen": max_tlen, 
            "count":count,
            "commonDesc":commonDesc,
            "rep_protein_id":rep_protein_id,
            "target_len" : int(lenMode)
            
        }
        newSameDomainList.append(domainInfo)
    newSameDomainList = sorted(newSameDomainList, key=lambda k: k["target_len"])
    architecture = {
        "totalSequence":totalSequence,
        "seqWithAnnotate":seqWithAnnotate,
        "architechure":newSameDomainList
    }
    return [tbl,architecture]

class HMMER:
                
    def parse_jackhmmer(self, hmm, domtbl, outFi, out_archiFi):
        align = parse_hmmer(hmm)
        tbl, architecture = parse_hmmerdomtbl(domtbl)
        data =[]
        for seq_id in tbl:
            info = tbl[seq_id]
            pw = align[seq_id]
            info['pairwise'] = pw
            data.append(info)
        
        dt_json = json.dumps(data)
        with open(outFi, 'w') as fo:
            json.dump(dt_json, fo)
        
        with open(out_archiFi, 'w') as fo_archi:
            json.dump( json.dumps(architecture), fo_archi)
            
class BLAST:
    def parse_psiblast(self, xml, outFi,out_archiFi):
        parse_psiblast_xml(xml,outFi,out_archiFi)


def parse_psiblast_xml(xml,outFi, out_archiFi):
    E_VALUE_THRESH = 0.001
    data = []
    sameDomain =defaultdict(list)
    seqWithAnnotate = 0

    with open(xml) as fh:
        records = NCBIXML.parse(fh)
        print(dir(records))
        count =0 
        for blast_record in records:
            query = blast_record.query
            query_sp = query.split()
            query_id = query_sp[0]
            count +=1
            print('---------------',count)
            for alignment in blast_record.alignments:
                num_matches = len(alignment.hsps)
                # print("num_matches:",num_matches)
                title = alignment.title
                title_split = title.split()
                title_id = title_split[0]
                title_seq_id = title_split[1]
                title_desc = title.split(']')[0] + ']'
                length = alignment.length
                # print(title_seq_id)
                align = ''
                for hsp in alignment.hsps:
                    desc = ''
                    if hsp.expect < E_VALUE_THRESH:
                        score = hsp.score
                        bits = hsp.bits
                        identities = hsp.identities
                        query_start = hsp.query_start
                        sbjct_start = hsp.sbjct_start
                        e_value = hsp.expect
                        query = hsp.query
                        match = hsp.match
                        sbjct = hsp.sbjct
                        align_length = hsp.align_length
                        # print(hsp)
                        # break
                        ident = round(identities / align_length * 100, 2)
                        gap_per = round(hsp.gaps / align_length * 100, 2)
                        query_cover = round(align_length / blast_record.query_length * 100, 2)
                       
                        align += 'hsp: ' + str(
                            alignment.hsps.index(hsp) + 1) + ': ' + str(hsp.sbjct_start) + ' to ' + str(
                            hsp.sbjct_start + align_length - 1) + ','
                            
                        align += f" Score = {str(score)}  bits ({str(hsp.bits)}), Expect = {str(e_value)} " 
                        align +=f"Identities = {identities}/{align_length} ({ident}%), Gaps = {hsp.gaps}/{align_length}  ({gap_per}%)\n\n"
                        #Positives = {hsp.gaps}/{align_length} (75%)
                        loops = round(align_length // 60 + 1)
                #        
                        for i in range(0, loops):
                            align += "Query %8s  %s %s\n" % (str(query_start), query[(i * 60):((i + 1) * 60)],
                                                             str(query_start + len(query[(i * 60):((i + 1) * 60)]) - 1))
                            align += "      %8s  %s\n" % ('', match[(i * 60):((i + 1) * 60)])
                            align += "Sbjct %8s  %s %s\n\n" % (str(sbjct_start), sbjct[(i * 60):((i + 1) * 60)], str(
                                sbjct_start + len(sbjct[(i * 60):((i + 1) * 60)]) - 1))
                            query_start += 60
                            sbjct_start += 60

                        cdd_names,cdd_ids, cdd_notes = retrieve_cdd(title_seq_id )
                        dt = {
                        "target":title_seq_id,
                        "query": query,
                        "query_covery": query_cover,
                        "target_len": align_length,
                        "E_value": e_value,
                        "score": float(score),
                        "domain": 1,
                        "domain_count": 1,
                        "acc": float(ident),
                        'desc': title_desc ,
                       "cdd_nameCat" :cdd_names ,
                        "cdd_noteCat" :cdd_notes,
                        "cdd_idCat" : cdd_ids,
                        "pairwise":[align]
                        }
                        data.append(dt)
                        target = title_seq_id
                        desc = title_desc
                        tlen = align_length
                        sameDomain[cdd_names].append([target, tlen, desc])
        print("number of round:",count)
        newSameDomainList = []
        totalSequence = len(data)
        for cdd_nameCat in sameDomain:
            min_tlen, max_tlen = 10000000,0
            descDt = defaultdict(int)
            proteinLenDt = defaultdict(list)
            count = len(sameDomain[cdd_nameCat])
            for target, tlen, desc in sameDomain[cdd_nameCat]:
                if cdd_nameCat != "":
                    seqWithAnnotate +=1
                proteinLenDt[str(tlen)].append(target)
                min_tlen = min(tlen,min_tlen)
                max_tlen =max(tlen,max_tlen)
                descDt[desc.split('[')[0].split(',')[0] ] +=1
            commonDesc = sorted(descDt.items(), key=lambda k: k[1], reverse=True)[0][0]
            lenMode = sorted(proteinLenDt.keys(), key=lambda k: int(k), reverse=True)[0]
            rep_protein_id = proteinLenDt[lenMode][0]
            domainInfo = {
                "cdd_nameCat": cdd_nameCat,
                "min_tlen": min_tlen,
                "max_tlen": max_tlen, 
                "count":count,
                "commonDesc":commonDesc,
                "rep_protein_id":rep_protein_id,
                "target_len" : int(lenMode)
                
            }
            newSameDomainList.append(domainInfo)
        newSameDomainList = sorted(newSameDomainList, key=lambda k: k["target_len"])
        architecture = {
            "totalSequence":totalSequence,
            "seqWithAnnotate":seqWithAnnotate,
            "architechure":newSameDomainList
        }
        with open(outFi, 'w') as fo, open(out_archiFi, 'w') as fo_archi:
            dt_json = json.dumps(data)
            json.dump(dt_json, fo)
            json.dump( json.dumps(architecture), fo_archi)

if __name__ == "__main__":
    fire.Fire(HMMER)
