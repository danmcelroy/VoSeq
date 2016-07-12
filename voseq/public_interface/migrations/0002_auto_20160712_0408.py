# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='date_extraction',
            field=models.DateField(null=True, blank=True),
        ),
    ]
