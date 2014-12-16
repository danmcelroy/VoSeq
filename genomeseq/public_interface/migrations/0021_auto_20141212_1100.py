# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0020_auto_20141210_1231'),
    ]

    operations = [
        migrations.AlterField(
            model_name='sequences',
            name='code',
            field=models.ForeignKey(to='public_interface.Vouchers'),
        ),
        migrations.AlterField(
            model_name='sequences',
            name='time_created',
            field=models.DateTimeField(null=True, auto_now_add=True),
        ),
        migrations.AlterField(
            model_name='sequences',
            name='time_edited',
            field=models.DateTimeField(null=True, auto_now=True),
        ),
    ]
