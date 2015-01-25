# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0014_auto_20150105_0106'),
    ]

    operations = [
        migrations.AlterField(
            model_name='sequences',
            name='time_created',
            field=models.DateTimeField(null=True, auto_now_add=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='sequences',
            name='time_edited',
            field=models.DateTimeField(null=True, auto_now=True),
            preserve_default=True,
        ),
    ]
