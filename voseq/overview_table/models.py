from django.db import models

from public_interface.models import Vouchers


class OverviewTable(Vouchers):
    sequences = models.TextField(help_text="HTML string of cells with length of"
                                           "sequences for each gene.")
