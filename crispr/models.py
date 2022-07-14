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

    def __str__(self):
        return self.gene_name + '_' + self.protein_name


class AlignTMScore(models.Model):
    chain1 = models.CharField(max_length=100, db_index=True)
    chain2_acc = models.CharField(max_length=100, db_index=True)
    chain1_len = models.IntegerField()
    chain2_len = models.IntegerField()
    align_len = models.IntegerField()
    cov_1 = models.IntegerField()
    cov_2 = models.IntegerField()
    RMSD = models.IntegerField()
    seq_ID = models.IntegerField()
    TMscore_1 = models.IntegerField()
    TMscore_2 = models.IntegerField()
    d0_1 = models.IntegerField()
    d0_2 = models.IntegerField()


class AlignSPScore(models.Model):
    chain1 = models.CharField(max_length=200, db_index=True)
    chain2_acc = models.CharField(max_length=200, db_index=True)
    tools = models.CharField(max_length=100)
    chain1_len = models.IntegerField()
    chain2_len = models.IntegerField()
    align_len = models.IntegerField()
    cov_1 = models.IntegerField()
    cov_2 = models.IntegerField()
    RMSD = models.IntegerField()
    seq_ID = models.IntegerField()
    TMscore_1 = models.IntegerField()
    TMscore_2 = models.IntegerField()
    d0_1 = models.IntegerField()
    d0_2 = models.IntegerField()
