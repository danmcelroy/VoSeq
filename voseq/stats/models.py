from django.db import models


class Stats(models.Model):
    """
    Statistics about number of vouchers, photos, etc to show in page footer.
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
