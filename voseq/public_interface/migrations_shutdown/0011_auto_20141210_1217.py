# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0010_auto_20141210_1215'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='auctor',
            field=models.CharField(max_length=100, null=True, blank=True, help_text='Person that described this taxon.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='country',
            field=models.CharField(max_length=100, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='determinedBy',
            field=models.CharField(max_length=100, null=True, blank=True, help_text='Person that identified the taxon for this specimen.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='edits',
            field=models.TextField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='family',
            field=models.CharField(max_length=100, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='hostorg',
            field=models.CharField(max_length=200, null=True, blank=True, help_text='Hostplant or other host.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='latesteditor',
            field=models.TextField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='notes',
            field=models.TextField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='orden',
            field=models.CharField(max_length=100, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='publishedIn',
            field=models.TextField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='species',
            field=models.CharField(max_length=100, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='subfamily',
            field=models.CharField(max_length=100, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='subtribe',
            field=models.CharField(max_length=100, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='tribe',
            field=models.CharField(max_length=100, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='voucherCode',
            field=models.CharField(max_length=100, null=True, blank=True, help_text='Original code of voucher specimen.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='voucherLocality',
            field=models.CharField(max_length=200, null=True, blank=True),
        ),
    ]
