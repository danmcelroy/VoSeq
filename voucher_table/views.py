from django.shortcuts import render
from django.http import HttpResponseRedirect

from .forms import VoucherTableForm
from .utils import VoucherTable
from core.utils import get_context


def index(request):
    context = get_context(request)
    context["form"] = VoucherTableForm()
    return render(request, 'voucher_table/index.html', context)


def results(request):
    context = get_context(request)
    if request.method == 'POST':
        form = VoucherTableForm(request.POST)
        if form.is_valid():
            table = VoucherTable(form.cleaned_data)
            response = table.create_csv_file()
            return response
        else:
            context["form"] = form
            return render(request, 'voucher_table/index.html', context)
    return HttpResponseRedirect('/create_voucher_table/')
