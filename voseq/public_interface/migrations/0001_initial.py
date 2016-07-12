# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='FlickrImages',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('voucher_image', models.URLField(help_text='URLs of the Flickr page.', blank=True)),
                ('thumbnail', models.URLField(help_text='URLs for the small sized image from Flickr.')),
                ('flickr_id', models.CharField(max_length=100, help_text='ID numbers from Flickr for our photo.')),
                ('image_file', models.ImageField(upload_to='', help_text='Placeholder for image file so we can send it to Flickr. The file has been deleted right after upload.', blank=True)),
            ],
            options={
                'verbose_name_plural': 'Flickr Images',
            },
        ),
        migrations.CreateModel(
            name='Genes',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('gene_code', models.CharField(max_length=100)),
                ('genetic_code', models.PositiveSmallIntegerField(help_text='Translation table (as number). See <a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi">http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi</a>', null=True)),
                ('length', models.PositiveSmallIntegerField(help_text='Number of base pairs', null=True)),
                ('description', models.CharField(max_length=255, help_text='Long gene name.', blank=True)),
                ('reading_frame', models.PositiveSmallIntegerField(help_text='Either 1, 2 or 3', null=True)),
                ('notes', models.TextField(blank=True)),
                ('aligned', models.CharField(choices=[('yes', 'yes'), ('no', 'no'), ('notset', 'notset')], max_length=6, default='notset')),
                ('intron', models.CharField(max_length=255, blank=True)),
                ('prot_code', models.CharField(choices=[('yes', 'yes'), ('no', 'no'), ('notset', 'notset')], max_length=6, default='notset')),
                ('gene_type', models.CharField(max_length=255, help_text='Nuclear, mitochondrial.', blank=True)),
                ('time_created', models.DateTimeField(auto_now_add=True)),
            ],
            options={
                'verbose_name_plural': 'Genes',
            },
        ),
        migrations.CreateModel(
            name='GeneSets',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('geneset_name', models.CharField(max_length=75)),
                ('geneset_creator', models.CharField(max_length=75)),
                ('geneset_description', models.CharField(max_length=140, blank=True)),
                ('geneset_list', models.TextField(help_text='As items separated by linebreak.')),
            ],
            options={
                'verbose_name_plural': 'Gene sets',
            },
        ),
        migrations.CreateModel(
            name='LocalImages',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('voucher_image', models.ImageField(upload_to='', help_text='voucher photo.', blank=True)),
            ],
            options={
                'verbose_name_plural': 'Local Images',
            },
        ),
        migrations.CreateModel(
            name='Primers',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('primer_f', models.CharField(max_length=100, blank=True)),
                ('primer_r', models.CharField(max_length=100, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='Sequences',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('gene_code', models.CharField(max_length=100)),
                ('sequences', models.TextField(blank=True)),
                ('accession', models.CharField(max_length=100, blank=True)),
                ('lab_person', models.CharField(max_length=100, blank=True)),
                ('time_created', models.DateTimeField(auto_now_add=True, null=True)),
                ('time_edited', models.DateTimeField(auto_now=True, null=True)),
                ('notes', models.TextField(blank=True)),
                ('genbank', models.NullBooleanField()),
                ('total_number_bp', models.IntegerField(null=True, blank=True)),
                ('number_ambiguous_bp', models.IntegerField(null=True, blank=True)),
            ],
            options={
                'verbose_name_plural': 'Sequences',
            },
        ),
        migrations.CreateModel(
            name='TaxonSets',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('taxonset_name', models.CharField(max_length=75)),
                ('taxonset_creator', models.CharField(max_length=75)),
                ('taxonset_description', models.CharField(max_length=140, blank=True)),
                ('taxonset_list', models.TextField(help_text='As items separated by linebreak.')),
            ],
            options={
                'verbose_name_plural': 'Taxon sets',
            },
        ),
        migrations.CreateModel(
            name='Vouchers',
            fields=[
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('code', models.CharField(serialize=False, max_length=300, primary_key=True, help_text='Voucher code.', unique=True)),
                ('orden', models.TextField(blank=True)),
                ('superfamily', models.TextField(blank=True)),
                ('family', models.TextField(blank=True)),
                ('subfamily', models.TextField(blank=True)),
                ('tribe', models.TextField(blank=True)),
                ('subtribe', models.TextField(blank=True)),
                ('genus', models.TextField(blank=True)),
                ('species', models.TextField(blank=True)),
                ('subspecies', models.TextField(blank=True)),
                ('country', models.TextField(blank=True)),
                ('specific_locality', models.TextField(help_text='Locality of origin for this specimen.', blank=True)),
                ('type_species', models.CharField(choices=[('unknown', 'unknown'), ('yes', 'yes'), ('not', 'not')], max_length=100, help_text='Is this a type species?')),
                ('latitude', models.FloatField(null=True, blank=True)),
                ('longitude', models.FloatField(null=True, blank=True)),
                ('max_altitude', models.IntegerField(help_text='Enter altitude in meters above sea level.', null=True, blank=True)),
                ('min_altitude', models.IntegerField(help_text='Enter altitude in meters above sea level.', null=True, blank=True)),
                ('collector', models.TextField(blank=True)),
                ('date_collection', models.DateField(null=True)),
                ('extraction', models.TextField(help_text='Number of extraction event.', blank=True)),
                ('extraction_tube', models.TextField(help_text='Tube containing DNA extract.', blank=True)),
                ('date_extraction', models.DateField(null=True)),
                ('extractor', models.TextField(blank=True)),
                ('voucher_locality', models.TextField(blank=True)),
                ('published_in', models.TextField(null=True, blank=True)),
                ('notes', models.TextField(null=True, blank=True)),
                ('edits', models.TextField(null=True, blank=True)),
                ('latest_editor', models.TextField(null=True, blank=True)),
                ('hostorg', models.TextField(help_text='Hostplant or other host.', blank=True)),
                ('sex', models.CharField(choices=[('male', 'male'), ('female', 'female'), ('larva', 'larva'), ('worker', 'worker'), ('queen', 'queen'), ('unknown', 'unknown')], max_length=100, blank=True)),
                ('voucher', models.CharField(choices=[('spread', 'spread'), ('in envelope', 'in envelope'), ('only photo', 'only photo'), ('no voucher', 'no voucher'), ('destroyed', 'destroyed'), ('lost', 'lost'), ('unknown', 'unknown')], max_length=100, help_text='Voucher status.', blank=True)),
                ('voucher_code', models.TextField(help_text='Alternative code of voucher specimen.', blank=True)),
                ('code_bold', models.TextField(help_text='Optional code for specimens kept in the BOLD database.', blank=True)),
                ('determined_by', models.TextField(help_text='Person that identified the taxon for this specimen.', blank=True)),
                ('author', models.TextField(help_text='Person that described this taxon.', blank=True)),
            ],
            options={
                'verbose_name_plural': 'Vouchers',
            },
        ),
        migrations.AddField(
            model_name='sequences',
            name='code',
            field=models.ForeignKey(to='public_interface.Vouchers', help_text='This is your voucher code.'),
        ),
        migrations.AddField(
            model_name='primers',
            name='for_sequence',
            field=models.ForeignKey(to='public_interface.Sequences', help_text='relation to Sequences table with reference for code and gene_code.'),
        ),
        migrations.AddField(
            model_name='localimages',
            name='voucher',
            field=models.ForeignKey(to='public_interface.Vouchers', help_text='Relation with id of voucher.'),
        ),
        migrations.AddField(
            model_name='flickrimages',
            name='voucher',
            field=models.ForeignKey(to='public_interface.Vouchers', help_text='Relation with id of voucher. Save as lower case.'),
        ),
    ]
