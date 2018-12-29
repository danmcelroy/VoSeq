import logging

from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import render

from core.utils import get_context
from public_interface.forms.admin_forms import VoucherForm


log = logging.getLogger(__name__)


@login_required
def add_voucher(request):
    context = get_context(request)
    form = VoucherForm(request.POST or None)
    if form.is_valid():
        voucher = form.save()
        voucher.user = request.user
        voucher.save()
        return HttpResponse("form is valid")

    context["form"] = form
    return render(request, 'admin_interface/add_voucher.html', context)
