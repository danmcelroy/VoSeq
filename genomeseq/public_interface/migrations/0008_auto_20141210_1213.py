# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0007_auto_20141210_1205'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='subspecies',
            field=models.CharField(max_length=100, blank=True),
        ),
    ]
