# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0019_auto_20150118_1340'),
    ]

    operations = [
        migrations.AlterField(
            model_name='genes',
            name='genetic_code',
            field=models.PositiveSmallIntegerField(blank=True, null=True, help_text='Translation table'),
            preserve_default=True,
        ),
    ]
