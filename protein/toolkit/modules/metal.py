from mymetal import mbp

class Mymetal:

    def __init__(self,fasta) -> None:
        self.fasta = fasta
    
    def run(self):
        a = mbp.predict(self.fasta)


class Main:
    def run(self,fasta):
        from mymetal import mbp,iof
        metal = Mymetal(fasta)
        iof.save_out_csv(metal, 'out_filename.csv')