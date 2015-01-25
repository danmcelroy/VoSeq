# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0003_auto_20141218_1352'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='typeSpecies',
            field=models.CharField(default='-', null=True, choices=[('-', '-'), ('d', "don't know"), ('y', 'yes'), ('n', 'no')], help_text='Is this a type species?', max_length=1, blank=True),
        ),
    ]
