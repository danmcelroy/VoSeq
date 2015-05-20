from django.shortcuts import render
from django.http import HttpResponseRedirect

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
    if request.method == 'POST':
        form = VoucherTableForm(request.POST)
        if form.is_valid():
            print(form.cleaned_data)
            table = VoucherTable(form.cleaned_data)
            response = table.create_csv_file()
            return response

    return HttpResponseRedirect('/create_voucher_table/')
