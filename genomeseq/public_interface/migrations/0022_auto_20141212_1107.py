# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0021_auto_20141212_1100'),
    ]

    operations = [
        migrations.AlterField(
            model_name='sequences',
            name='genbank',
            field=models.NullBooleanField(),
        ),
        migrations.AlterField(
            model_name='sequences',
            name='sequences',
            field=models.TextField(blank=True),
        ),
    ]
