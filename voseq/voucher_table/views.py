from django.shortcuts import render
from django.http import HttpResponseRedirect

from .forms import VoucherTableForm
from .utils import VoucherTable
from core.utils import get_version_stats
from core.utils import get_username


def index(request):
    VERSION, STATS = get_version_stats()
    return render(
        request,
        'voucher_table/index.html',
        {
            'username': get_username(request),
            'version': VERSION,
            'stats': STATS,
            'form': VoucherTableForm(),
        },
    )


def results(request):
    VERSION, STATS = get_version_stats()
    if request.method == 'POST':
        form = VoucherTableForm(request.POST)
        if form.is_valid():
            table = VoucherTable(form.cleaned_data)
            response = table.create_csv_file()
            return response
        else:
            return render(
                request,
                'voucher_table/index.html',
                {
                    'username': get_username(request),
                    'version': VERSION,
                    'stats': STATS,
                    'form': form,
                }
            )
    return HttpResponseRedirect('/create_voucher_table/')
