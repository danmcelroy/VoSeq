# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Stats',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('vouchers', models.IntegerField(help_text='Number of records, or vouchers.')),
                ('orders', models.IntegerField(help_text='Number of Orders.')),
                ('families', models.IntegerField()),
                ('genera', models.IntegerField()),
                ('species', models.IntegerField()),
                ('sequences', models.IntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='VouchersPerGene',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('gene_code', models.CharField(max_length=100)),
                ('voucher_count', models.IntegerField()),
            ],
        ),
    ]
