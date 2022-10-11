from django.db import models

# Create your models here.

class CASInfo(models.Model):

    cas_class = models.CharField(max_length=100, db_index=True)
    accession = models.CharField(max_length=100, db_index=True)
    entry_name = models.CharField(max_length=100, db_index=True)
    data_class = models.CharField(max_length=100, db_index=True)
    protein_name = models.CharField(max_length=1000, db_index=True)
    gene_name = models.CharField(max_length=200, db_index=True)
    organism = models.CharField(max_length=1000, db_index=True)
    taxonomy_id = models.CharField(max_length=100, null=True)
    sequence_length = models.IntegerField(null=True)
    genome_genbank  = models.CharField(max_length=200, null=True)
    protein_genebankID = models.CharField(max_length=200, null=True)

    def __str__(self):
        return self.gene_name + '_' + self.protein_name


class AlignTMScore(models.Model):
    chain1 = models.CharField(max_length=100, db_index=True)
    chain2_acc = models.CharField(max_length=100, db_index=True)
    chain1_len = models.IntegerField()
    chain2_len = models.IntegerField()
    align_len = models.IntegerField()
    cov_1 = models.FloatField()
    cov_2 = models.FloatField()
    RMSD = models.FloatField()
    seq_ID = models.FloatField()
    TMscore_1 = models.FloatField()
    TMscore_2 = models.FloatField()
    d0_1 = models.FloatField()
    d0_2 = models.FloatField()


class AlignSPScore(models.Model):
    chain1 = models.CharField(max_length=200, db_index=True)
    chain2_acc = models.CharField(max_length=200, db_index=True)
    chain1_len = models.IntegerField()
    chain2_len = models.IntegerField()
    tool = models.CharField(max_length=100)
    eff_len = models.IntegerField()
    RMSD = models.FloatField()
    SPe = models.FloatField()
    SPa = models.FloatField()
    SPb = models.FloatField()
    seq_ID = models.FloatField(null=True)


class AlignFatcatScore(models.Model):
    chain1 = models.CharField(max_length=200, db_index=True)
    chain2_acc = models.CharField(max_length=200, db_index=True)
    chain1_len = models.IntegerField()
    chain2_len = models.IntegerField()
    cov1 = models.FloatField()
    cov2 = models.FloatField()
    seq_ID = models.FloatField()
    similar = models.FloatField()
    alignScore = models.FloatField()
    tmScore = models.FloatField()
    RMSD = models.FloatField()
