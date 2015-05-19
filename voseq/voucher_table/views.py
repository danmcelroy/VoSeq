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
from create_dataset.utils import CreateDataset
from public_interface.models import Genes


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
            table = GeneTable(form.cleaned_data)
            response = create_excel_file(table.stats)
            return response

    return render(request, 'gene_table/index.html',
                  {
                      'version': version,
                      'stats': stats,
                      'form': VoucherTableForm(),
                  },
                  )
