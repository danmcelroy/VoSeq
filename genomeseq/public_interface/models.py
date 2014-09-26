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
