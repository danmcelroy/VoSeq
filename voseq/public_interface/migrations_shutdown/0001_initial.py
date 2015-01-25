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
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('gene_code', models.CharField(max_length=255)),
                ('genetic_code', models.PositiveSmallIntegerField(blank=True)),
                ('length', models.PositiveSmallIntegerField()),
                ('description', models.CharField(max_length=255, blank=True)),
                ('reading_frame', models.PositiveSmallIntegerField()),
                ('notes', models.TextField(blank=True)),
                ('aligned', models.CharField(choices=[('yes', 'yes'), ('no', 'no'), ('notset', 'notset')], max_length=6, default='notset')),
                ('intron', models.CharField(max_length=255, blank=True)),
                ('prot_code', models.CharField(choices=[('yes', 'yes'), ('no', 'no'), ('notset', 'notset')], max_length=6, default='notset')),
                ('gene_type', models.CharField(max_length=255, blank=True)),
                ('time_created', models.DateTimeField(auto_now_add=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='GeneSets',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('geneset_name', models.CharField(max_length=75, default=None)),
                ('geneset_creator', models.CharField(max_length=75, default=None)),
                ('geneset_description', models.CharField(max_length=100, blank=True, default=None)),
                ('geneset_list', models.TextField(default=None)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Members',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('firstname', models.CharField(max_length=100)),
                ('lastname', models.CharField(max_length=100)),
                ('login', models.CharField(max_length=100)),
                ('passwd', models.CharField(max_length=100)),
                ('admin', models.BinaryField(default=None)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Primers',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('code', models.CharField(max_length=100)),
                ('gene_code', models.CharField(max_length=100)),
                ('primer1', models.CharField(max_length=100, blank=True)),
                ('primer2', models.CharField(max_length=100, blank=True)),
                ('primer3', models.CharField(max_length=100, blank=True)),
                ('primer4', models.CharField(max_length=100, blank=True)),
                ('primer5', models.CharField(max_length=100, blank=True)),
                ('primer6', models.CharField(max_length=100, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Sequences',
            fields=[
                ('id', models.AutoField(verbose_name='ID', auto_created=True, primary_key=True, serialize=False)),
                ('code', models.CharField(max_length=100)),
                ('gene_code', models.CharField(max_length=100)),
                ('sequences', models.TextField()),
                ('accession', models.CharField(max_length=100, blank=True)),
                ('labPerson', models.CharField(max_length=100, blank=True)),
                ('time_created', models.DateTimeField(auto_now_add=True)),
                ('time_edited', models.DateTimeField(auto_now=True)),
                ('notes', models.TextField(blank=True)),
                ('genbank', models.BooleanField(default=None)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
