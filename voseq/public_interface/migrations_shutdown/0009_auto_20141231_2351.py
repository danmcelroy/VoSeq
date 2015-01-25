# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0008_auto_20141231_2348'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='max_altitude',
            field=models.IntegerField(blank=True, null=True, help_text='Enter altitude in meters above sea level.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='min_altitude',
            field=models.IntegerField(blank=True, null=True, help_text='Enter altitude in meters above sea level.'),
        ),
    ]
