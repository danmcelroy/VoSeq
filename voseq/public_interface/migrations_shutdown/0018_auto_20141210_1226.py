# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0017_auto_20141210_1223'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='extraction',
            field=models.CharField(blank=True, max_length=50, null=True, help_text='Number of extraction event.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='extractionTube',
            field=models.CharField(blank=True, max_length=50, null=True, help_text='Tube containing DNA extract.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='timestamp',
            field=models.DateTimeField(null=True),
        ),
    ]
