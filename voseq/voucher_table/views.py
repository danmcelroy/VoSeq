from collections import OrderedDict
import csv
import os
import uuid

from django.shortcuts import render
from django.http import HttpResponse
from amas import AMAS

from .forms import VoucherTableForm
from core.utils import get_version_stats
from core.utils import get_voucher_codes
from core.utils import get_gene_codes
from public_interface.models import Vouchers
from create_dataset.utils import CreateDataset


def index(request):
    version, stats = get_version_stats()
    form = VoucherTableForm()

    return render(request, 'voucher_table/index.html',
                  {
                      'version': version,
                      'stats': stats,
                      'form': form,
                  },
                  )


def results(request):
    version, stats = get_version_stats()
    if request.method == 'POST':
        form = VoucherTableForm(request.POST)
        if form.is_valid():
            table = VoucherTable(form.cleaned_data)
            # response = create_excel_file(table.stats)
            return response

    return render(request, 'gene_table/index.html',
                  {
                      'version': version,
                      'stats': stats,
                      'form': VoucherTableForm(),
                  },
                  )


class VoucherTable(object):
    def __init__(self, cleaned_data):
        print(cleaned_data)
        self.cleaned_data = cleaned_data
        self.voucher_codes = get_voucher_codes(cleaned_data)
        self.gene_codes = get_gene_codes(cleaned_data)
        self.voucher_info_values = self.get_voucher_info_values()

    def get_voucher_info_values(self):
        return [i.lower() for i in self.cleaned_data['voucher_info']]

    def get_voucher_info(self):
        vouchers = Vouchers.objects.all().values(*self.voucher_info_values)
        print(">>> cleaned_data", self.cleaned_data)
        print(">>> vouchers", vouchers)
