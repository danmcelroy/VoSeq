# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0009_auto_20141210_1214'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='max_altitude',
            field=models.IntegerField(blank=True, help_text='Enter altitude in meters above sea level.', null=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='min_altitude',
            field=models.IntegerField(blank=True, help_text='Enter altitude in meters above sea level.', null=True),
        ),
    ]
