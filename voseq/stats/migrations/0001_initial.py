# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Stats',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, verbose_name='ID', serialize=False)),
                ('vouchers', models.IntegerField(help_text='Number of records, or vouchers.')),
                ('orders', models.IntegerField(help_text='Number of Orders.')),
                ('families', models.IntegerField()),
                ('genera', models.IntegerField()),
                ('species', models.IntegerField()),
                ('sequences', models.IntegerField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
