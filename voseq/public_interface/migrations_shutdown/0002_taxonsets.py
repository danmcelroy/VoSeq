# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='TaxonSets',
            fields=[
                ('id', models.AutoField(auto_created=True, serialize=False, verbose_name='ID', primary_key=True)),
                ('taxonset_name', models.CharField(max_length=50)),
                ('taxonset_creator', models.CharField(max_length=75)),
                ('taxonset_description', models.CharField(max_length=100)),
                ('taxonset_list', models.TextField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
