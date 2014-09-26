# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Genes',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('gene_code', models.CharField(max_length=255)),
                ('genetic_code', models.PositiveSmallIntegerField(blank=True)),
                ('length', models.PositiveSmallIntegerField()),
                ('description', models.CharField(blank=True, max_length=255)),
                ('reading_frame', models.PositiveSmallIntegerField()),
                ('notes', models.TextField(blank=True)),
                ('aligned', models.CharField(choices=[('yes', 'yes'), ('no', 'no'), ('notset', 'notset')], max_length=6, default='notset')),
                ('intron', models.CharField(blank=True, max_length=255)),
                ('prot_code', models.CharField(choices=[('yes', 'yes'), ('no', 'no'), ('notset', 'notset')], max_length=6, default='notset')),
                ('gene_type', models.CharField(blank=True, max_length=255)),
                ('time_created', models.DateTimeField(auto_now_add=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
