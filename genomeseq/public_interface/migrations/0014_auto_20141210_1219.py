# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0013_auto_20141210_1219'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='extractor',
            field=models.CharField(max_length=100, blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='specificLocality',
            field=models.CharField(help_text='Locality of origin for this specimen.', max_length=100, blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='typeSpecies',
            field=models.CharField(help_text='Is this a type species?', max_length=1, choices=[('d', 'd'), ('y', 'y'), ('n', 'n')], blank=True, null=True),
        ),
    ]
