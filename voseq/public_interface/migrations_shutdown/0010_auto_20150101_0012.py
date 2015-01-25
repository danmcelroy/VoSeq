# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0009_auto_20141231_2351'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='collector',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='country',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='extractionTube',
            field=models.CharField(max_length=50, help_text='Tube containing DNA extract.', blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='family',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='genus',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='orden',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='species',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='subfamily',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='subtribe',
            field=models.CharField(max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='tribe',
            field=models.CharField(max_length=100, blank=True),
        ),
    ]
