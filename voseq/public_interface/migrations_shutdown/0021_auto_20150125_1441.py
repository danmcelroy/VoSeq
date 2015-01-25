# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0020_auto_20150119_2231'),
    ]

    operations = [
        migrations.CreateModel(
            name='FlickrImage',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('voucherImage', models.URLField(help_text='URLs of the Flickr page.')),
                ('thumbnail', models.URLField(help_text='URLs for the small sized image from Flickr.')),
                ('flickr_id', models.CharField(help_text='ID numbers from Flickr for our photo.', max_length=100)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Primer',
            fields=[
                ('id', models.AutoField(auto_created=True, verbose_name='ID', serialize=False, primary_key=True)),
                ('primer_f', models.CharField(max_length=100, blank=True)),
                ('primer_r', models.CharField(max_length=100, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Voucher',
            fields=[
                ('code', models.CharField(help_text='Voucher code.', unique=True, max_length=100, primary_key=True, serialize=False)),
                ('orden', models.CharField(max_length=100, blank=True)),
                ('superfamily', models.CharField(max_length=100, blank=True)),
                ('family', models.CharField(max_length=100, blank=True)),
                ('subfamily', models.CharField(max_length=100, blank=True)),
                ('tribe', models.CharField(max_length=100, blank=True)),
                ('subtribe', models.CharField(max_length=100, blank=True)),
                ('genus', models.CharField(max_length=100, blank=True)),
                ('species', models.CharField(max_length=100, blank=True)),
                ('subspecies', models.CharField(max_length=100, blank=True)),
                ('country', models.CharField(max_length=100, blank=True)),
                ('specificLocality', models.CharField(help_text='Locality of origin for this specimen.', max_length=250, blank=True)),
                ('typeSpecies', models.CharField(choices=[('d', "don't know"), ('y', 'yes'), ('n', 'no')], help_text='Is this a type species?', max_length=1)),
                ('latitude', models.FloatField(null=True, blank=True)),
                ('longitude', models.FloatField(null=True, blank=True)),
                ('max_altitude', models.IntegerField(help_text='Enter altitude in meters above sea level.', null=True, blank=True)),
                ('min_altitude', models.IntegerField(help_text='Enter altitude in meters above sea level.', null=True, blank=True)),
                ('collector', models.CharField(max_length=100, blank=True)),
                ('dateCollection', models.DateField(null=True)),
                ('extraction', models.CharField(help_text='Number of extraction event.', max_length=50, blank=True)),
                ('extractionTube', models.CharField(help_text='Tube containing DNA extract.', max_length=50, blank=True)),
                ('dateExtraction', models.DateField(null=True)),
                ('extractor', models.CharField(max_length=100, blank=True)),
                ('voucherLocality', models.CharField(max_length=200, blank=True)),
                ('publishedIn', models.TextField(null=True, blank=True)),
                ('notes', models.TextField(null=True, blank=True)),
                ('edits', models.TextField(null=True, blank=True)),
                ('latesteditor', models.TextField(null=True, blank=True)),
                ('hostorg', models.CharField(help_text='Hostplant or other host.', max_length=200, blank=True)),
                ('sex', models.CharField(choices=[('m', 'male'), ('f', 'female'), ('l', 'larva'), ('w', 'worker'), ('q', 'queen')], max_length=1, blank=True)),
                ('voucher', models.CharField(choices=[('s', 'spread'), ('e', 'in envelope'), ('p', 'only photo'), ('n', 'no voucher'), ('d', 'destroyed'), ('l', 'lost')], max_length=1, blank=True)),
                ('voucherCode', models.CharField(help_text='Original code of voucher specimen.', max_length=100, blank=True)),
                ('code_bold', models.CharField(help_text='Optional code for specimens kept in the BOLD database.', max_length=100, blank=True)),
                ('determinedBy', models.CharField(help_text='Person that identified the taxon for this specimen.', max_length=100, blank=True)),
                ('auctor', models.CharField(help_text='Person that described this taxon.', max_length=100, blank=True)),
                ('timestamp', models.DateTimeField(null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.RenameModel(
            old_name='Genes',
            new_name='Gene',
        ),
        migrations.RenameModel(
            old_name='GeneSets',
            new_name='GeneSet',
        ),
        migrations.RenameModel(
            old_name='Members',
            new_name='Member',
        ),
        migrations.RenameModel(
            old_name='Sequences',
            new_name='Sequence',
        ),
        migrations.RenameModel(
            old_name='TaxonSets',
            new_name='TaxonSet',
        ),
        migrations.RemoveField(
            model_name='flickrimages',
            name='voucher',
        ),
        migrations.DeleteModel(
            name='FlickrImages',
        ),
        migrations.RemoveField(
            model_name='primers',
            name='for_sequence',
        ),
        migrations.DeleteModel(
            name='Primers',
        ),
        migrations.DeleteModel(
            name='Vouchers',
        ),
        migrations.AddField(
            model_name='primer',
            name='for_sequence',
            field=models.ForeignKey(to='public_interface.Sequence', help_text='relation to Sequence table with reference for code and gene_code.'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='flickrimage',
            name='voucher',
            field=models.ForeignKey(to='public_interface.Voucher', help_text='Relation with id of voucher. Save as lower case.'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='sequence',
            name='code',
            field=models.ForeignKey(to='public_interface.Voucher', help_text='Save as lower case.'),
            preserve_default=True,
        ),
    ]
