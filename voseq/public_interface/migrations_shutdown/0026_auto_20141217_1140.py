# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0025_auto_20141217_1136'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='typeSpecies',
            field=models.CharField(help_text='Is this a type species?', blank=True, null=True, choices=[('d', "don't know"), ('y', 'yes'), ('n', 'no')], max_length=1),
        ),
    ]
