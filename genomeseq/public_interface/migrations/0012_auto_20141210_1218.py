# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0011_auto_20141210_1217'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='sex',
            field=models.CharField(choices=[('m', 'male'), ('f', 'female'), ('l', 'larva'), ('w', 'worker'), ('q', 'queen')], max_length=1, blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='voucher',
            field=models.CharField(choices=[('s', 'spread'), ('e', 'in envelope'), ('p', 'only photo'), ('n', 'no voucher'), ('d', 'destroyed'), ('l', 'lost')], max_length=1, blank=True, null=True),
        ),
    ]
