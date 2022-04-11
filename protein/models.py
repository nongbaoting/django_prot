from unittest import result
from django.db import models
from django.utils import timezone
# Create your models here.
import datetime

# from django_prot import protein

class SubmitInfo(models.Model):
    email_addr = models.EmailField(max_length=50)
    proj_name = models.CharField(max_length=50)
    job_name = models.CharField(max_length=50, null=True)
    tools = models.CharField(max_length=50)
    params = models.TextField(max_length=10000, null=True)
    upload_date = models.DateTimeField('date published')
    running_date = models.DateTimeField('date published', null=True)
    completed_date = models.DateTimeField('date published', null=True)
    run_status = models.CharField(max_length=50, null=True)

    def __str__(self):
        line = self.job_name
        if self.run_status is not None:
            line += "\t" + self.run_status
        return line


class GeneInfo(models.Model):
    entry_name = models.CharField(max_length=1000, null=True)
    gene_name = models.CharField(max_length=2550, null=True)
    accession = models.CharField(max_length=255, null=True)
    status = models.CharField(max_length=255, null=True)


class Structure(models.Model):
    accession = models.ForeignKey(GeneInfo, on_delete=models.CASCADE)
    alphfold = models.CharField(max_length=255)


class AF_uniprot(models.Model):
    uniprot = models.CharField(max_length=255)
    uniprot_name = models.CharField(max_length=2550)
    status = models.CharField(max_length=255)
    protein_names = models.CharField(max_length=2550)
    gene_names = models.CharField(max_length=2550)
    organism = models.CharField(max_length=255)
    length = models.IntegerField()

# blast
class SubmitInfoNew(models.Model):
    uuid = models.CharField(max_length=50)
    proj_type = models.CharField(max_length=50)
    job_name = models.CharField(max_length=50, null=True)
    email_addr = models.EmailField(max_length=50, null=True)
    tools = models.CharField(max_length=1000)
    params = models.TextField(max_length=10000, null=True)
    upload_date = models.DateTimeField('date published')
    running_date = models.DateTimeField('date published', null=True)
    completed_date = models.DateTimeField('date published', null=True)
    run_status = models.CharField(max_length=255, null=True)
    task_status = models.CharField(max_length=20, null=True)
   
    def __str__(self):
        line = self.job_name 
        if self.run_status is not None:
            line += "\t" + self.run_status
        return line

class CDD(models.Model):
    pssm = models.CharField(max_length=100, db_index=True,)
    cdd_id = models.CharField(max_length=100, db_index=True, unique=True, primary_key=True)
    cdd_name = models.CharField(max_length=255, db_index=True)
    cdd_desc =  models.TextField()
    cdd_desc_short = models.CharField(max_length=1000,db_index=True, default='')
    pssm_len = models.IntegerField(db_index=True,)  

class NrInfo(models.Model):
    protin_id = models.CharField(max_length=200, db_index=True)
    desc  = models.TextField()
    length = models.IntegerField(db_index=True,)
    filename = models.CharField(max_length=200, db_index=True)

    def __str__(self):
        return self.protin_id + "\t" + self.desc

class protCD(models.Model):
    cdd_id = models.ForeignKey(CDD, related_name= "cdd_protCD", on_delete=models.CASCADE)
    cdd_name = models.CharField(max_length=255, db_index=True)
    protin_id = models.CharField(max_length=200, db_index=True)
    length = models.IntegerField(db_index=True,)  
    pssm  = models.CharField(max_length=100, db_index=True,)
    start = models.IntegerField(db_index=True,)  
    end = models.IntegerField(db_index=True,)  
    evalue = models.CharField(max_length=100, db_index=True,)
    biscore = models.CharField(max_length=100, db_index=True,)
    
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_name

class protCDncbi(models.Model):
    cdd_id = models.ForeignKey(CDD, related_name= "cdd_protCDncbi", on_delete=models.CASCADE)
    cdd_name = models.CharField(max_length=255, db_index=True)
    protin_id =  models.CharField(max_length=200, db_index=True)
    length = models.IntegerField(db_index=True,)  
    pssm  = models.CharField(max_length=100, db_index=True,)
    start = models.IntegerField(db_index=True,)  
    end = models.IntegerField(db_index=True,)  
    evalue = models.CharField(max_length=100, db_index=True,)
    biscore = models.CharField(max_length=100, db_index=True,)
   
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_name

class protCDncbiCat(models.Model):
    protin_id =  models.ForeignKey(NrInfo, related_name= "nrinfo_protCDncbi", on_delete=models.CASCADE)
    cdd_nameCat = models.TextField( )
    cdd_idCat = models.TextField(  )
    
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_nameCat

class protCDncbiOne(models.Model):
    protin_id = models.CharField(max_length=200, db_index=True)
    protin_nr =  models.ForeignKey(NrInfo, related_name= "nrinfo_protCDncbiOne", on_delete=models.CASCADE)
    cdd_nameCat = models.TextField( )
    cdd_idCat = models.TextField(  )
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_nameCat

class NrCD(models.Model):
    cdd_id = models.ForeignKey(CDD, related_name= "cdd_nrcd", on_delete=models.CASCADE)
    protin_id =  models.CharField(max_length=200, db_index=True)
    length = models.IntegerField(db_index=True,)  
    desc  = models.CharField(max_length=1000, db_index=True,)
    pssm  = models.CharField(max_length=100, db_index=True,)
    start = models.IntegerField(db_index=True,)  
    end = models.IntegerField(db_index=True,)  
    evalue = models.CharField(max_length=100, db_index=True,)
    biscore = models.CharField(max_length=100, db_index=True,)
    cdd_name = models.CharField(max_length=1000,db_index=True,)
    cdd_nameCat = models.TextField( )
    cdd_idCat = models.TextField( null=True )
    
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_name

class ProtCDD(models.Model):
    protin_id =  models.CharField(max_length=200, db_index=True, unique=True, primary_key=True)
    length = models.IntegerField(db_index=True,)
    cdd_ids = models.TextField()
    cdd_annots = models.TextField()
    cdd_notes = models.TextField( )
    cdd_locs = models.TextField( )
    
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_annots

class ProtCDDeach(models.Model):
    protin_id =  models.CharField(max_length=200, db_index=True, )
    length = models.IntegerField(db_index=True,)
    cdd_id = models.CharField(max_length=200,db_index=True,)
    cdd_name = models.CharField(max_length=1000,db_index=True,)
    cdd_note = models.CharField(max_length=1000, db_index=True,)
    cdd_loca = models.CharField(max_length=200,db_index=True,) 
    
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_name

class ProtCDDeach2(models.Model):
    protin_id =  models.CharField(max_length=200, db_index=True, )
    length = models.IntegerField(db_index=True,)
    cdd_id = models.CharField(max_length=200,db_index=True,)
    cdd_name = models.CharField(max_length=1000,db_index=True,)
    cdd_note = models.CharField(max_length=1000, db_index=True,)
    cdd_loca = models.CharField(max_length=200,db_index=True,) 
    
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_name
    
class ProtCDD2_head(models.Model):
    protin_id =  models.CharField(max_length=200, db_index=True, unique=True, primary_key=True)
    length = models.IntegerField(db_index=True,)
    cdd_ids = models.TextField()
    cdd_annots = models.TextField()
    cdd_notes = models.TextField( )
    cdd_locs = models.TextField( )
    
    def __str__(self):
        return self.protin_id + "\t" + self.cdd_annots

