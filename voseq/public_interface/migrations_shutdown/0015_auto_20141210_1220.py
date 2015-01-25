# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0014_auto_20141210_1219'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='collector',
            field=models.CharField(blank=True, null=True, max_length=100),
        ),
    ]
