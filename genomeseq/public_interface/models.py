from django.db import models


class Genes(models.Model):
    gene_code = models.CharField(max_length=255)
    genetic_code = models.PositiveSmallIntegerField(blank=True)
    length = models.PositiveSmallIntegerField()
    description = models.CharField(max_length=255, blank=True)
    reading_frame = models.PositiveSmallIntegerField()
    notes = models.TextField(blank=True)
    aligned = models.CharField(max_length=6, choices=(
        ('yes', 'yes'),
        ('no', 'no'),
        ('notset', 'notset'),
    ),
                               default='notset')
    intron = models.CharField(max_length=255, blank=True)
    prot_code = models.CharField(max_length=6, choices=(
        ('yes', 'yes'),
        ('no', 'no'),
        ('notset', 'notset'),
    ),
                               default='notset')
    gene_type = models.CharField(max_length=255, blank=True)

    time_created = models.DateTimeField(auto_now_add=True)


class GeneSets(models.Model):
    geneset_name = models.CharField(max_length=75, default=None)
    geneset_creator = models.CharField(max_length=75, default=None)
    geneset_description = models.CharField(max_length=100, default=None, blank=True)
    geneset_list = models.TextField(default=None)


class Members(models.Model):
    firstname = models.CharField(max_length=100)
    lastname = models.CharField(max_length=100)
    login = models.CharField(max_length=100)
    passwd = models.CharField(max_length=100)
    admin = models.BinaryField(default=None)


class Primers(models.Model):
    code = models.CharField(max_length=100)
    gene_code = models.CharField(max_length=100)
    primer1 = models.CharField(max_length=100, blank=True)
    primer2 = models.CharField(max_length=100, blank=True)
    primer3 = models.CharField(max_length=100, blank=True)
    primer4 = models.CharField(max_length=100, blank=True)
    primer5 = models.CharField(max_length=100, blank=True)
    primer6 = models.CharField(max_length=100, blank=True)


class Sequences(models.Model):
    code = models.CharField(max_length=100)
    gene_code = models.CharField(max_length=100)
    sequences = models.TextField()
    accession = models.CharField(max_length=100, blank=True)
    labPerson = models.CharField(max_length=100, blank=True)
    time_created = models.DateTimeField(auto_now_add=True)
    time_edited = models.DateTimeField(auto_now=True)
    notes = models.TextField(blank=True)
    genbank = models.BooleanField(default=None, blank=True)


class TaxonSets(models.Model):
    taxonset_name = models.CharField(max_length=50)
    taxonset_creator = models.CharField(max_length=75)
    taxonset_description = models.CharField(max_length=100)
    taxonset_list = models.TextField()



