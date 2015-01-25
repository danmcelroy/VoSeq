# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0011_auto_20150101_0752'),
    ]

    operations = [
        migrations.AlterField(
            model_name='flickrimages',
            name='voucher',
            field=models.ForeignKey(help_text='Relation with id of voucher. Save as lower case.', to='public_interface.Vouchers'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='determinedBy',
            field=models.CharField(max_length=100, help_text='Person that identified the taxon for this specimen.', blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='extraction',
            field=models.CharField(max_length=50, help_text='Number of extraction event.', blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='extractor',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='sex',
            field=models.CharField(max_length=1, blank=True, choices=[('m', 'male'), ('f', 'female'), ('l', 'larva'), ('w', 'worker'), ('q', 'queen')]),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='voucher',
            field=models.CharField(max_length=1, blank=True, choices=[('s', 'spread'), ('e', 'in envelope'), ('p', 'only photo'), ('n', 'no voucher'), ('d', 'destroyed'), ('l', 'lost')]),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='voucherCode',
            field=models.CharField(max_length=100, help_text='Original code of voucher specimen.', blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='voucherLocality',
            field=models.CharField(max_length=200, blank=True),
        ),
    ]
