# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0001_squashed_0026_auto_20141217_1140'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='superfamily',
            field=models.CharField(max_length=100, null=True),
        ),
    ]
