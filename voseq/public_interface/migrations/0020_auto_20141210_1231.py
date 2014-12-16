# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0019_auto_20141210_1228'),
    ]

    operations = [
        migrations.AlterField(
            model_name='flickrimages',
            name='voucher',
            field=models.ForeignKey(to='public_interface.Vouchers', help_text='Relation with id of voucher', null=True, blank=True),
        ),
    ]
