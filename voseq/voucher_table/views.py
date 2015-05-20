from django.shortcuts import render

from .forms import VoucherTableForm
from .utils import VoucherTable
from core.utils import get_version_stats


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
