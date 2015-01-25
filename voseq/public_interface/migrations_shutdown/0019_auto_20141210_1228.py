# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0018_auto_20141210_1226'),
    ]

    operations = [
        migrations.AlterField(
            model_name='flickrimages',
            name='flickr_id',
            field=models.CharField(max_length=100, help_text='ID numbers from Flickr for our photo.'),
        ),
    ]
