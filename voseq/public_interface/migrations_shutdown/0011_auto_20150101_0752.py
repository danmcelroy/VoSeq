# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0010_auto_20150101_0012'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='specificLocality',
            field=models.CharField(max_length=250, help_text='Locality of origin for this specimen.', blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='typeSpecies',
            field=models.CharField(max_length=1, help_text='Is this a type species?', choices=[('d', "don't know"), ('y', 'yes'), ('n', 'no')]),
        ),
    ]
