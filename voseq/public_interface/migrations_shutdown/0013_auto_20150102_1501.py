# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0012_auto_20150102_0600'),
    ]

    operations = [
        migrations.RenameField(
            model_name='primers',
            old_name='primer1',
            new_name='primer_f',
        ),
        migrations.RenameField(
            model_name='primers',
            old_name='primer2',
            new_name='primer_r',
        ),
        migrations.RemoveField(
            model_name='primers',
            name='code',
        ),
        migrations.RemoveField(
            model_name='primers',
            name='gene_code',
        ),
        migrations.RemoveField(
            model_name='primers',
            name='primer3',
        ),
        migrations.RemoveField(
            model_name='primers',
            name='primer4',
        ),
        migrations.RemoveField(
            model_name='primers',
            name='primer5',
        ),
        migrations.RemoveField(
            model_name='primers',
            name='primer6',
        ),
        migrations.AddField(
            model_name='primers',
            name='for_sequence',
            field=models.ForeignKey(help_text='relation to Sequences table with reference for code and gene_code.', to='public_interface.Sequences', default=''),
            preserve_default=False,
        ),
    ]
