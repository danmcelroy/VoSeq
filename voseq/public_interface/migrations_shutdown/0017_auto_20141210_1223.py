# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0016_auto_20141210_1222'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='genus',
            field=models.CharField(max_length=100, null=True, blank=True),
        ),
    ]
