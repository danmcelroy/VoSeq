# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.utils.timezone import utc
import datetime


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0013_auto_20150102_1501'),
    ]

    operations = [
        migrations.AlterField(
            model_name='sequences',
            name='time_created',
            field=models.DateTimeField(default=datetime.datetime(2015, 1, 5, 1, 6, 14, 266202, tzinfo=utc), editable=False),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='sequences',
            name='time_edited',
            field=models.DateTimeField(default=datetime.datetime(2015, 1, 5, 1, 6, 27, 553956, tzinfo=utc)),
            preserve_default=False,
        ),
    ]
