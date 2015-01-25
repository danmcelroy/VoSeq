# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0008_auto_20141210_1213'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='subspecies',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
    ]
