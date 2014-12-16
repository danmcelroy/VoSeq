# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0015_auto_20141210_1220'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='specificLocality',
            field=models.CharField(max_length=250, blank=True, null=True, help_text='Locality of origin for this specimen.'),
        ),
    ]
