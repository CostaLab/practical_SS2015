from django.db import models

# Create your models here.
class Gene(models.Model):
    gene_name = models.CharField(max_length=200)
    refFile = models.FileField(upload_to='documents',default='------')
    
    def __str__(self):
        return self.gene_name

class Technology(models.Model):
    gene = models.ForeignKey(Gene,related_name='technologies')
    technology = models.CharField(max_length=200) 
    sequence = models.CharField(max_length=500,default='000000')
    motif_file = models.FileField(upload_to='documents',default='------')
    snp_file = models.FileField(upload_to='documents', null=True, blank=True)
    barGraph = models.ImageField(upload_to='graphics',default='------')
    compared_with = models.CharField(max_length=500,default='000000')
    commonGraph = models.ImageField(upload_to='graphics',default='------')
    differGraph = models.ImageField(upload_to='graphics', default='-----')
    
    def __str__(self):
        return self.technology