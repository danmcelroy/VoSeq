# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0022_auto_20141212_1107'),
    ]

    operations = [
        migrations.AddField(
            model_name='vouchers',
            name='superfamily',
            field=models.CharField(blank=True, max_length=100, null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='flickrimages',
            name='voucher',
            field=models.ForeignKey(help_text='Relation with id of voucher. Save as lower case.', to='public_interface.Vouchers', blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='sequences',
            name='code',
            field=models.ForeignKey(to='public_interface.Vouchers', help_text='Save as lower case.'),
        ),
    ]
