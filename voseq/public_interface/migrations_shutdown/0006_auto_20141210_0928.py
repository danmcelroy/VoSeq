# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0005_auto_20141210_0917'),
    ]

    operations = [
        migrations.CreateModel(
            name='FlickrImages',
            fields=[
                ('id', models.AutoField(serialize=False, primary_key=True, verbose_name='ID', auto_created=True)),
                ('voucherImage', models.URLField(help_text='URLs of the Flickr page.')),
                ('thumbnail', models.URLField(help_text='URLs for the small sized image from Flickr.')),
                ('flickr_id', models.IntegerField(help_text='ID numbers from Flickr for our photo.')),
                ('voucher', models.ForeignKey(to='public_interface.Vouchers', help_text='Relation with id of voucher')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.RemoveField(
            model_name='vouchers',
            name='flickr_id',
        ),
        migrations.RemoveField(
            model_name='vouchers',
            name='thumbnail',
        ),
        migrations.RemoveField(
            model_name='vouchers',
            name='voucherImage',
        ),
    ]
