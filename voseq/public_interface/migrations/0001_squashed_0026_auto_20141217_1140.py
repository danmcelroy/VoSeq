# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    replaces = [('public_interface', '0001_initial'), ('public_interface', '0002_taxonsets'), ('public_interface', '0003_vouchers'), ('public_interface', '0004_auto_20141210_0835'), ('public_interface', '0005_auto_20141210_0917'), ('public_interface', '0006_auto_20141210_0928'), ('public_interface', '0007_auto_20141210_1205'), ('public_interface', '0008_auto_20141210_1213'), ('public_interface', '0009_auto_20141210_1214'), ('public_interface', '0010_auto_20141210_1215'), ('public_interface', '0011_auto_20141210_1217'), ('public_interface', '0012_auto_20141210_1218'), ('public_interface', '0013_auto_20141210_1219'), ('public_interface', '0014_auto_20141210_1219'), ('public_interface', '0015_auto_20141210_1220'), ('public_interface', '0016_auto_20141210_1222'), ('public_interface', '0017_auto_20141210_1223'), ('public_interface', '0018_auto_20141210_1226'), ('public_interface', '0019_auto_20141210_1228'), ('public_interface', '0020_auto_20141210_1231'), ('public_interface', '0021_auto_20141212_1100'), ('public_interface', '0022_auto_20141212_1107'), ('public_interface', '0023_auto_20141212_1433'), ('public_interface', '0024_auto_20141216_1157'), ('public_interface', '0025_auto_20141217_1136'), ('public_interface', '0026_auto_20141217_1140')]

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Vouchers',
            fields=[
                ('code', models.CharField(unique=True, primary_key=True, max_length=100, help_text='Voucher code.')),
                ('orden', models.CharField(max_length=100)),
                ('superfamily', models.CharField(max_length=100, blank=True)),
                ('family', models.CharField(max_length=100)),
                ('subfamily', models.CharField(max_length=100)),
                ('tribe', models.CharField(max_length=100)),
                ('subtribe', models.CharField(max_length=100)),
                ('genus', models.CharField(max_length=100)),
                ('species', models.CharField(max_length=100)),
                ('subspecies', models.CharField(max_length=100)),
                ('country', models.CharField(max_length=100)),
                ('specificLocality', models.CharField(max_length=100, help_text='Locality of origin for this specimen.')),
                ('typeSpecies', models.CharField(max_length=100, help_text='Is this a type species?')),
                ('latitude', models.FloatField()),
                ('longitude', models.FloatField()),
                ('max_altitude', models.IntegerField(help_text='Enter altitude in meters above sea level.')),
                ('collector', models.CharField(max_length=100)),
                ('dateCollection', models.DateField(null=True)),
                ('voucherImage', models.TextField(help_text='URLs of the Flickr page.')),
                ('thumbnail', models.TextField(help_text='URLs for the small sized image from Flickr.')),
                ('extraction', models.CharField(max_length=50, help_text='Number of extraction event.')),
                ('extractionTube', models.CharField(max_length=50, help_text='Tube containing DNA extract.')),
                ('dateExtraction', models.DateField(null=True)),
                ('extractor', models.CharField(max_length=100)),
                ('voucherLocality', models.CharField(max_length=200)),
                ('publishedIn', models.TextField()),
                ('notes', models.TextField()),
                ('edits', models.TextField()),
                ('latesteditor', models.TextField()),
                ('hostorg', models.CharField(max_length=200, help_text='Hostplant or other host.')),
                ('sex', models.CharField(choices=[('m', 'male'), ('f', 'female'), ('l', 'larva'), ('w', 'worker'), ('q', 'queen')], max_length=1)),
                ('voucher', models.CharField(choices=[('s', 'spread'), ('e', 'in envelope'), ('p', 'only photo'), ('n', 'no voucher')], max_length=1)),
                ('voucherCode', models.CharField(max_length=100, help_text='Original code of voucher specimen.')),
                ('flickr_id', models.TextField(help_text='ID numbers from Flickr for our photo.')),
                ('determinedBy', models.CharField(max_length=100, help_text='Person that identified the taxon for this specimen.')),
                ('auctor', models.CharField(max_length=100, help_text='Person that described this taxon.')),
                ('timestamp', models.DateTimeField()),
                ('min_altitude', models.IntegerField(help_text='Enter altitude in meters above sea level.', default=None)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Genes',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('gene_code', models.CharField(max_length=255)),
                ('genetic_code', models.PositiveSmallIntegerField(blank=True)),
                ('length', models.PositiveSmallIntegerField()),
                ('description', models.CharField(blank=True, max_length=255)),
                ('reading_frame', models.PositiveSmallIntegerField()),
                ('notes', models.TextField(blank=True)),
                ('aligned', models.CharField(choices=[('yes', 'yes'), ('no', 'no'), ('notset', 'notset')], max_length=6, default='notset')),
                ('intron', models.CharField(blank=True, max_length=255)),
                ('prot_code', models.CharField(choices=[('yes', 'yes'), ('no', 'no'), ('notset', 'notset')], max_length=6, default='notset')),
                ('gene_type', models.CharField(blank=True, max_length=255)),
                ('time_created', models.DateTimeField(auto_now_add=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='GeneSets',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('geneset_name', models.CharField(max_length=75, default=None)),
                ('geneset_creator', models.CharField(max_length=75, default=None)),
                ('geneset_description', models.CharField(blank=True, max_length=100, default=None)),
                ('geneset_list', models.TextField(default=None)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Members',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('firstname', models.CharField(max_length=100)),
                ('lastname', models.CharField(max_length=100)),
                ('login', models.CharField(max_length=100)),
                ('passwd', models.CharField(max_length=100)),
                ('admin', models.BinaryField(default=None)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Primers',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('code', models.CharField(max_length=100)),
                ('gene_code', models.CharField(max_length=100)),
                ('primer1', models.CharField(blank=True, max_length=100)),
                ('primer2', models.CharField(blank=True, max_length=100)),
                ('primer3', models.CharField(blank=True, max_length=100)),
                ('primer4', models.CharField(blank=True, max_length=100)),
                ('primer5', models.CharField(blank=True, max_length=100)),
                ('primer6', models.CharField(blank=True, max_length=100)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Sequences',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('code', models.ForeignKey(to='public_interface.Vouchers', help_text='Save as lower case.')),
                ('gene_code', models.CharField(max_length=100)),
                ('sequences', models.TextField(blank=True)),
                ('accession', models.CharField(blank=True, max_length=100)),
                ('labPerson', models.CharField(blank=True, max_length=100)),
                ('time_created', models.DateTimeField(auto_now_add=True, null=True)),
                ('time_edited', models.DateTimeField(auto_now=True, null=True)),
                ('notes', models.TextField(blank=True)),
                ('genbank', models.NullBooleanField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TaxonSets',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('taxonset_name', models.CharField(max_length=50)),
                ('taxonset_creator', models.CharField(max_length=75)),
                ('taxonset_description', models.CharField(max_length=100)),
                ('taxonset_list', models.TextField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='FlickrImages',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, verbose_name='ID', primary_key=True)),
                ('voucherImage', models.URLField(help_text='URLs of the Flickr page.')),
                ('thumbnail', models.URLField(help_text='URLs for the small sized image from Flickr.')),
                ('flickr_id', models.CharField(max_length=100, help_text='ID numbers from Flickr for our photo.')),
                ('voucher', models.ForeignKey(to='public_interface.Vouchers', null=True, blank=True, help_text='Relation with id of voucher. Save as lower case.')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.RemoveField(
            model_name='Vouchers',
            name='flickr_id',
        ),
        migrations.RemoveField(
            model_name='Vouchers',
            name='thumbnail',
        ),
        migrations.RemoveField(
            model_name='Vouchers',
            name='voucherImage',
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='typeSpecies',
            field=models.CharField(choices=[('d', 'd'), ('y', 'y'), ('n', 'n')], max_length=1, help_text='Is this a type species?'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='voucher',
            field=models.CharField(choices=[('s', 'spread'), ('e', 'in envelope'), ('p', 'only photo'), ('n', 'no voucher'), ('d', 'destroyed'), ('l', 'lost')], max_length=1),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='subspecies',
            field=models.CharField(blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='subspecies',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='max_altitude',
            field=models.IntegerField(null=True, blank=True, help_text='Enter altitude in meters above sea level.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='min_altitude',
            field=models.IntegerField(null=True, blank=True, help_text='Enter altitude in meters above sea level.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='auctor',
            field=models.CharField(null=True, blank=True, max_length=100, help_text='Person that described this taxon.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='country',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='determinedBy',
            field=models.CharField(null=True, blank=True, max_length=100, help_text='Person that identified the taxon for this specimen.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='edits',
            field=models.TextField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='family',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='hostorg',
            field=models.CharField(null=True, blank=True, max_length=200, help_text='Hostplant or other host.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='latesteditor',
            field=models.TextField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='notes',
            field=models.TextField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='orden',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='publishedIn',
            field=models.TextField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='species',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='subfamily',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='subtribe',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='tribe',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='voucherCode',
            field=models.CharField(null=True, blank=True, max_length=100, help_text='Original code of voucher specimen.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='voucherLocality',
            field=models.CharField(null=True, blank=True, max_length=200),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='sex',
            field=models.CharField(choices=[('m', 'male'), ('f', 'female'), ('l', 'larva'), ('w', 'worker'), ('q', 'queen')], null=True, blank=True, max_length=1),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='voucher',
            field=models.CharField(choices=[('s', 'spread'), ('e', 'in envelope'), ('p', 'only photo'), ('n', 'no voucher'), ('d', 'destroyed'), ('l', 'lost')], null=True, blank=True, max_length=1),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='latitude',
            field=models.FloatField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='longitude',
            field=models.FloatField(null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='extractor',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='specificLocality',
            field=models.CharField(null=True, blank=True, max_length=100, help_text='Locality of origin for this specimen.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='typeSpecies',
            field=models.CharField(choices=[('d', 'd'), ('y', 'y'), ('n', 'n')], null=True, blank=True, max_length=1, help_text='Is this a type species?'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='collector',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='specificLocality',
            field=models.CharField(null=True, blank=True, max_length=250, help_text='Locality of origin for this specimen.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='genus',
            field=models.CharField(null=True, blank=True, max_length=100),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='extraction',
            field=models.CharField(null=True, blank=True, max_length=50, help_text='Number of extraction event.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='extractionTube',
            field=models.CharField(null=True, blank=True, max_length=50, help_text='Tube containing DNA extract.'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='timestamp',
            field=models.DateTimeField(null=True),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='typeSpecies',
            field=models.CharField(choices=[('d', "don't know"), ('y', 'yes'), ('n', 'no')], null=True, blank=True, max_length=1, help_text='Is this a type species?'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='typeSpecies',
            field=models.CharField(choices=[('d', "don't know"), ('y', 'yes'), ('n', 'no')], default='', null=True, blank=True, max_length=1, help_text='Is this a type species?'),
        ),
        migrations.AlterField(
            model_name='Vouchers',
            name='typeSpecies',
            field=models.CharField(choices=[('d', "don't know"), ('y', 'yes'), ('n', 'no')], null=True, blank=True, max_length=1, help_text='Is this a type species?'),
        ),
    ]
