# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0002_taxonsets'),
    ]

    operations = [
        migrations.CreateModel(
            name='Vouchers',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True, serialize=False)),
                ('code', models.CharField(unique=True, help_text='Voucher code.', max_length=100)),
                ('orden', models.CharField(max_length=100)),
                ('family', models.CharField(max_length=100)),
                ('subfamily', models.CharField(max_length=100)),
                ('tribe', models.CharField(max_length=100)),
                ('subtribe', models.CharField(max_length=100)),
                ('genus', models.CharField(max_length=100)),
                ('species', models.CharField(max_length=100)),
                ('subspecies', models.CharField(max_length=100)),
                ('country', models.CharField(max_length=100)),
                ('specificLocality', models.CharField(help_text='Locality of origin for this specimen.', max_length=100)),
                ('typeSpecies', models.CharField(help_text='Is this a type species?', max_length=100)),
                ('latitude', models.FloatField()),
                ('longitude', models.FloatField()),
                ('altitude', models.IntegerField(help_text='Enter altitude in meters above sea level.')),
                ('collector', models.CharField(max_length=100)),
                ('dateCollection', models.DateField(null=True)),
                ('voucherImage', models.URLField(help_text='URL of the Flickr page.')),
                ('thumbnail', models.URLField(help_text='URL for the small sized image from Flickr.')),
                ('extraction', models.CharField(help_text='Number of extraction event.', max_length=50)),
                ('extractionTube', models.CharField(help_text='Tube containing DNA extract.', max_length=50)),
                ('dateExtraction', models.DateField(null=True)),
                ('extractor', models.CharField(max_length=100)),
                ('voucherLocality', models.CharField(max_length=200)),
                ('publishedIn', models.TextField()),
                ('notes', models.TextField()),
                ('edits', models.TextField()),
                ('latesteditor', models.TextField()),
                ('hostorg', models.CharField(help_text='Hostplant or other host.', max_length=200)),
                ('sex', models.CharField(max_length=1, choices=[('m', 'male'), ('f', 'female'), ('l', 'larva'), ('w', 'worker'), ('q', 'queen')])),
                ('voucher', models.CharField(max_length=1, choices=[('s', 'spread'), ('e', 'in envelope'), ('p', 'only photo'), ('n', 'no voucher')])),
                ('voucherCode', models.CharField(help_text='Original code of voucher specimen.', max_length=100)),
                ('flickr_id', models.IntegerField(help_text='ID number from Flickr for our photo.')),
                ('determinedBy', models.CharField(help_text='Person that identified the taxon for this specimen.', max_length=100)),
                ('auctor', models.CharField(help_text='Person that described this taxon.', max_length=100)),
                ('timestamp', models.DateTimeField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
