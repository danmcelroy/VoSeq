# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0006_auto_20141210_0928'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='vouchers',
            name='id',
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='code',
            field=models.CharField(serialize=False, primary_key=True, help_text='Voucher code.', max_length=100, unique=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='typeSpecies',
            field=models.CharField(choices=[('d', 'd'), ('y', 'y'), ('n', 'n')], max_length=1, help_text='Is this a type species?'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='voucher',
            field=models.CharField(choices=[('s', 'spread'), ('e', 'in envelope'), ('p', 'only photo'), ('n', 'no voucher'), ('d', 'destroyed'), ('l', 'lost')], max_length=1),
        ),
    ]
