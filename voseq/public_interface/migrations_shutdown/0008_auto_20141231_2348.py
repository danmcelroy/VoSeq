# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0007_auto_20141231_2345'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='auctor',
            field=models.CharField(help_text='Person that described this taxon.', max_length=100, blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='max_altitude',
            field=models.IntegerField(help_text='Enter altitude in meters above sea level.', blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='min_altitude',
            field=models.IntegerField(help_text='Enter altitude in meters above sea level.', blank=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='superfamily',
            field=models.CharField(max_length=100, blank=True),
        ),
    ]
