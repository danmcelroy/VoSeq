# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0015_auto_20150105_0114'),
    ]

    operations = [
        migrations.AddField(
            model_name='vouchers',
            name='code_bold',
            field=models.CharField(max_length=100, help_text='Optional code for specimens kept in the BOLD database.', blank=True),
            preserve_default=True,
        ),
    ]
