
import re
import os
from collections import defaultdict
from crispr.utils.myFunctions import scanAndFind_pattern, calAlignIdent, run_cmd
from Bio import PDB
from protein.toolkit import myFunctions


class SPalign:

    def __init__(self, path_1, path_2, tempDir, outFi):
        self.path_1 = path_1
        self.path_2 = path_2
        self.tempDir = tempDir
        self.matrix = []
        self.logFi = ''
        self.outFi = outFi

    def align(self):
        self.logFi = self.outFi + "_Spalgin.txt"
        cmd = 'cd %s;SP-align.intel %s %s >%s' % (
            self.tempDir, self.path_1, self.path_2, self.logFi)
        run_cmd(cmd)

    def parse_LogFile(self):
        dt = defaultdict(dict)
        # reg_log = re.compile(
        #     r'RMSD/Nali= (?P<RMSD>\d+(\.\d+)?) / (?P<Nali>\d+) ;GDT= (?P<GDT>\d+(\.\d+)?) ;TMscore\(a,b,c\)= (?P<TMscore_a>\d+(\.\d+)?) (?P<TMscore_b>\d+(\.\d+)?) (?P<TMscore_c>\d+(\.\d+)?) SEQID= (?P<SEQID>\d+(\.\d+)?)%',
        #     re.MULTILINE)
        reg_log = re.compile(r'RMSD/Nali= (?P<RMSD>\d+(\.\d+)?) / (?P<Nali>\d+) ;GDT= (?P<GDT>\d+(\.\d+)?) ;TMscore\(a,b,c\)= (?P<TMscore_a>\d+(\.\d+)?) (?P<TMscore_b>\d+(\.\d+)?) (?P<TMscore_c>\d+(\.\d+)?)',
                             re.MULTILINE)
        # reg_so = re.compile(r"Structure Overlap: (?P<so>\d+(\.\d+)?) %")
        reg_sp = re.compile(
            r'################## Alignment Report ##################')
        reg_position = re.compile(
            r"\(':' denotes the residue pairs of distance <= 4A, and '.' denotes <=8A\)\n(?P<source_start>\d+)\s+(?P<seq_1>.*)\s+(?P<source_end>\d+)\n(?P<pairwise>.*)\n(?P<target_start>\d+)\s+(?P<seq_2>.*)\s+(?P<target_end>\d+)", re.M)
        reg_pFold = re.compile(
            r'Pfold= (?P<Pfold>\d+(?:\.\d+)?) \%; SPe/SPa/SPb= (?P<SPe>\d+(?:\.\d+)?) (?P<SPa>\d+(?:\.\d+)?) (?P<SPb>\d+(?:\.\d+)?) ;Effective_Length: (?P<eff_len>\d+)')
        reg_len = re.compile(r': Length= (?P<pdb1_len>\d+) (?P<pdb2_len>\d+)')

        reports_log = open(self.logFi).read().strip()
        for report in reg_sp.split(reports_log)[1:]:
            # print(report)
            lines = report.split('\n')
            pdb1, pdb2 = [i.split('.pdb')[0]
                          for i in lines[1].split(':')[0].split(' ')]
            matchs = reg_log.search(report)

            m_pos = reg_position.search(report)
            seq_1, seq_2, pairwise = m_pos.group('seq_1'), m_pos.group(
                'seq_2'),  m_pos.group('pairwise')

            Pfold, SPe, SPa, SPb, eff_len = reg_pFold.search(report).groups()
            pdb1_len, pdb2_len = reg_len.search(report).groups()
            align4_len, ali_ident = calAlignIdent(
                seq_1, seq_2, pairwise[5:])
            eff_ident = round(align4_len / float(eff_len), 2)
            RMSD = matchs.group('RMSD')

            struc_overlap = ''

            self.results = [pdb1, pdb2, pdb1_len, pdb2_len,
                            SPe, SPa, SPb, eff_len,
                            matchs.group('RMSD'), matchs.group(
                                'Nali'), matchs.group('GDT'),
                            matchs.group('TMscore_a'), matchs.group(
                                'TMscore_b'),
                            matchs.group('TMscore_c'),
                            struc_overlap, str(eff_ident), str(ali_ident),
                            m_pos.group('source_start'), m_pos.group('source_end'), m_pos.group(
                                'target_start'), m_pos.group('target_end'),
                            m_pos.group('seq_1'), m_pos.group(
                                'seq_2'),  m_pos.group('pairwise')
                            ]

            self.content = {
                "input_pdb": pdb1,
                "target_pdb": pdb2,
                "len_1": pdb1_len,
                "len_2": pdb2_len,
                "SPe": SPe,
                "SPa": SPa,
                "SPb": SPe,
                "eff_len": eff_len,
                'ali_ident': ali_ident,
                "eff_ident": eff_ident,
                "seq_1": seq_1, "seq_2": seq_2, "pairwise": pairwise,
                "RMSD": RMSD
            }

    def getMatrix(self):
        re_ma = re.compile(r"Rotation Matrix:")
        re_sp = re.compile(r"\s+")
        matrix = []
        with open(self.logFi, 'r') as f:
            m = 0
            for li in f:
                if m == 1:
                    cell = re_sp.split(li.strip())
                    self.matrix.append([float(i) for i in cell[0:]])
                if re_ma.match(li):
                    m += 1
                if len(self.matrix) >= 3:
                    print(self.matrix)
                    break

    def rotate(self, ):

        self.content['out1_pdb'] = self.outFi
        parser = PDB.PDBParser()
        name = os.path.basename(self.path_1).split('.')[0]
        struct = parser.get_structure(name, self.path_1)
        for model in struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        x, y, z = atom.coord
                        X = self.matrix[0][0] + self.matrix[0][1] * x + \
                            self.matrix[0][2] * y + self.matrix[0][3] * z

                        Y = self.matrix[1][0] + self.matrix[1][1] * x + \
                            self.matrix[1][2] * y + self.matrix[1][3] * z

                        Z = self.matrix[2][0] + self.matrix[2][1] * x + \
                            self.matrix[2][2] * y + self.matrix[2][3] * z
                        atom.coord = [X, Y, Z]

        io = PDB.PDBIO()
        io.set_structure(struct)
        io.save(self.outFi)

    def save_pickle(self,):
        item_pickle = self.outFi + '.pickle'
        myFunctions.pickle_dumpObj2file(self.content, item_pickle)


def pair_align(path1, path2, tempDir, outFi_name):
    outFi = os.path.join(tempDir, outFi_name)
    spalign = SPalign(path1, path2, tempDir, outFi)
    spalign.align()
    spalign.parse_LogFile()
    spalign.getMatrix()
    spalign.rotate()
    spalign.save_pickle()

    return spalign
