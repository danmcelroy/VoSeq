# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0003_vouchers'),
    ]

    operations = [
        migrations.RenameField(
            model_name='vouchers',
            old_name='altitude',
            new_name='max_altitude',
        ),
        migrations.AddField(
            model_name='vouchers',
            name='min_altitude',
            field=models.IntegerField(help_text='Enter altitude in meters above sea level.', default=None),
            preserve_default=False,
        ),
    ]
