# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0018_auto_20150117_2110'),
    ]

    operations = [
        migrations.AlterField(
            model_name='genesets',
            name='geneset_description',
            field=models.CharField(max_length=140, blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='taxonsets',
            name='taxonset_description',
            field=models.CharField(max_length=140, blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='taxonsets',
            name='taxonset_name',
            field=models.CharField(max_length=75),
            preserve_default=True,
        ),
    ]
