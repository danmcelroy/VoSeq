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
from public_interface.models import Sequences
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
            response = table.create_csv_file()
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
        self.cleaned_data = cleaned_data
        self.voucher_codes = get_voucher_codes(cleaned_data)
        self.gene_codes = get_gene_codes(cleaned_data)
        self.voucher_info_values = self.get_voucher_info_values()
        self.voucher_info = self.get_voucher_info()
        self.sequences_info = self.get_sequence_info()
        self.warnings = []

    def get_voucher_info_values(self):
        voucher_info_values = self.cleaned_data['voucher_info'] + self.cleaned_data['collector_info']
        return voucher_info_values

    def get_voucher_info(self):
        vouchers_info = Vouchers.objects.all().values(*self.voucher_info_values)
        return self.convert_voucher_info_to_dict(vouchers_info)

    def convert_voucher_info_to_dict(self, vouchers_info):
        my_dict = OrderedDict()
        for voucher in vouchers_info:
            code = voucher['code']
            if code in self.voucher_codes:
                my_dict[code] = voucher
        return my_dict

    def get_sequence_info(self):
        seq_values = OrderedDict()
        seq_info = Sequences.objects.all().values('code', 'gene_code', 'sequences')

        for seq in seq_info:
            code = seq['code']
            gene_code = seq['gene_code']
            if code in self.voucher_codes and gene_code in self.gene_codes:
                if code not in seq_values:
                    seq_values[code] = OrderedDict()
                seq_values[code][gene_code] = len(seq['sequences'])
        return seq_values

    def create_csv_file(self):
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="gene_table.csv"'

        writer = csv.writer(response)
        row = self.get_headers()
        writer.writerow(row)

        for voucher_code in self.voucher_codes:
            row = [voucher_code]

            try:
                self.voucher_info[voucher_code]
            except KeyError:
                warning = 'We don\'t have voucher {} in our database.'.format(voucher_code)
                self.warnings.append(warning)
                continue

            item = self.voucher_info[voucher_code]
            for i in self.voucher_info_values:
                if i == 'code':
                    continue
                row.append(item[i])

            for gene_code in self.gene_codes:
                try:
                    row.append(self.sequences_info[voucher_code][gene_code])
                except KeyError:
                    warning = "We don't have sequences for {} and {}".format(gene_code, voucher_code)
                    self.warnings.append(warning)
                    row.append('')

            writer.writerow(row)
        return response

    def get_headers(self):
        row = tuple()

        for i in self.voucher_info_values:
            if i == 'specificLocality':
                row += ('Specific Locality',)
                continue
            row += (i.capitalize(),)

        for i in self.gene_codes:
            row += (i,)
        return row
