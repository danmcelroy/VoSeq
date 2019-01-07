# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='OverviewTable',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('sequence_string', models.TextField(help_text='HTML string of cells with length of sequences for each gene.')),
                ('o_code', models.CharField(max_length=300)),
                ('orden', models.TextField(blank=True)),
                ('superfamily', models.TextField(blank=True)),
                ('family', models.TextField(blank=True)),
                ('subfamily', models.TextField(blank=True)),
                ('genus', models.TextField(blank=True)),
                ('species', models.TextField(blank=True)),
            ],
        ),
    ]
