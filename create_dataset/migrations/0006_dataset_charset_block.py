# Generated by Django 2.2.13 on 2020-11-28 18:09

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('create_dataset', '0005_dataset_sister_dataset_id'),
    ]

    operations = [
        migrations.AddField(
            model_name='dataset',
            name='charset_block',
            field=models.TextField(null=True),
        ),
    ]