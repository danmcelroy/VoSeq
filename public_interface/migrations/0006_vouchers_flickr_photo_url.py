# Generated by Django 2.2.13 on 2021-03-14 21:47

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0005_auto_20210313_1919'),
    ]

    operations = [
        migrations.AddField(
            model_name='vouchers',
            name='flickr_photo_url',
            field=models.TextField(blank=True, null=True),
        ),
    ]