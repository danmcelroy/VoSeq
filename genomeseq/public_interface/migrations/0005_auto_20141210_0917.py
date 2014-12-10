# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0004_auto_20141210_0835'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vouchers',
            name='flickr_id',
            field=models.TextField(help_text='ID numbers from Flickr for our photo.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='thumbnail',
            field=models.TextField(help_text='URLs for the small sized image from Flickr.'),
        ),
        migrations.AlterField(
            model_name='vouchers',
            name='voucherImage',
            field=models.TextField(help_text='URLs of the Flickr page.'),
        ),
    ]
