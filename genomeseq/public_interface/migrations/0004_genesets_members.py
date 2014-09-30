# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0003_auto_20140930_0926'),
    ]

    operations = [
        migrations.CreateModel(
            name='GeneSets',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
                ('geneset_name', models.CharField(default=None, max_length=75)),
                ('geneset_creator', models.CharField(default=None, max_length=75)),
                ('geneset_description', models.CharField(default=None, max_length=100, blank=True)),
                ('geneset_list', models.TextField(default=None)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Members',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, auto_created=True, verbose_name='ID')),
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
    ]
