# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0006_auto_20141231_2322'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='hostorg',
            field=models.CharField(max_length=200, help_text='Hostplant or other host.', blank=True),
        ),
    ]
