import datetime

from haystack import indexes

from .models import Vouchers


class VoucherIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    code = indexes.CharField(model_attr='voucherCode')

    def get_model(self):
        return Vouchers

    # TODO change to time_edited, time_created with auto in tables and migrate_db script
    """
    def index_queryset(self, using=None):
        # Used when the entire index for model is updated.
        return self.get_model().objects.filter(time_created__lte=datetime.datetime.now())
    """
