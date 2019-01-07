from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('public_interface', '0004_auto_20181231_0208'),
    ]

    operations = [
        migrations.RunSQL(
            "UPDATE public_interface_genes SET gene_type=1 "
            "WHERE public_interface_genes.gene_type ILIKE 'mitochondrial';"
        ),
        migrations.RunSQL(
            "UPDATE public_interface_genes SET gene_type=2 "
            "WHERE public_interface_genes.gene_type ILIKE 'nuclear';"
        ),
        migrations.RunSQL(
            "UPDATE public_interface_genes SET gene_type=3 "
            "WHERE public_interface_genes.gene_type ILIKE 'ribosomal';"
        ),
    ]
