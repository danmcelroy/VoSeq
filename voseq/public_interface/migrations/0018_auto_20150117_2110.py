# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0017_auto_20150115_2235'),
    ]

    operations = [
        migrations.AlterField(
            model_name='genesets',
            name='geneset_creator',
            field=models.CharField(max_length=75),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='genesets',
            name='geneset_description',
            field=models.CharField(blank=True, max_length=100),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='genesets',
            name='geneset_list',
            field=models.TextField(),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='genesets',
            name='geneset_name',
            field=models.CharField(max_length=75),
            preserve_default=True,
        ),
    ]
