# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0012_auto_20141210_1218'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='latitude',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='longitude',
            field=models.FloatField(blank=True, null=True),
        ),
    ]
