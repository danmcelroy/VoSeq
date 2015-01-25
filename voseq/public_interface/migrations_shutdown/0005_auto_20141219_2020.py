# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0004_auto_20141218_1353'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='typeSpecies',
            field=models.CharField(max_length=1, help_text='Is this a type species?', blank=True, null=True, choices=[('d', "don't know"), ('y', 'yes'), ('n', 'no')]),
        ),
    ]
