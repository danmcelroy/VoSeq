# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0005_auto_20141219_2020'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='subspecies',
            field=models.CharField(blank=True, max_length=100),
        ),
    ]
