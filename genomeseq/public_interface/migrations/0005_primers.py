# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0004_genesets_members'),
    ]

    operations = [
        migrations.CreateModel(
            name='Primers',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('code', models.CharField(max_length=100)),
                ('gene_code', models.CharField(max_length=100)),
                ('primer1', models.CharField(blank=True, max_length=100)),
                ('primer2', models.CharField(blank=True, max_length=100)),
                ('primer3', models.CharField(blank=True, max_length=100)),
                ('primer4', models.CharField(blank=True, max_length=100)),
                ('primer5', models.CharField(blank=True, max_length=100)),
                ('primer6', models.CharField(blank=True, max_length=100)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
