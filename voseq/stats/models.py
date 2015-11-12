from django.db import models


class Stats(models.Model):
    """Statistics about number of vouchers, photos, etc to show in page footer.
    """
    vouchers = models.IntegerField(
        help_text='Number of records, or vouchers.'
    )
    orders = models.IntegerField(
        help_text='Number of Orders.'
    )
    families = models.IntegerField()
    genera = models.IntegerField()
    species = models.IntegerField()
    sequences = models.IntegerField()

    def __str__(self):
        return "{0} vouchers".format(self.vouchers)


class VouchersPerGene(models.Model):
    """Number of vouchers that have sequences for each of our genes.
    """
    gene_code = models.CharField(max_length=100)
    voucher_count = models.IntegerField()

    def __str__(self):
        return "{0} vouchers in gene {1}".format(self.voucher_count, self.gene_code)
