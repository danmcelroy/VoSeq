from datetime import datetime
import logging
from typing import Dict, Any

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpRequest
from django.shortcuts import render, redirect, reverse

from core.utils import get_context
from public_interface.forms.admin_forms import VoucherForm
from public_interface.models import Vouchers


log = logging.getLogger(__name__)


@login_required
def add_voucher(request: HttpRequest, code: str = None) -> HttpResponse:
    context = get_context(request)
    form = None
    voucher = None
    if code:
        voucher = Vouchers.objects.filter(code=code)
        if voucher:
            form = VoucherForm(request.POST or None, initial=voucher.values()[0])

    if not form:
        form = VoucherForm(request.POST or None)

    if form.is_valid():
        voucher = form.save()
        voucher.user = request.user
        voucher.save()
        return HttpResponse("form is valid")
    elif "update" in request.POST:
        # we are editing this voucher
        data = clean_post_data(request)
        voucher = Vouchers.objects.filter(code=request.POST.get("code"))
        voucher.update(**data)
        messages.add_message(request, messages.INFO,
                             f"The voucher of code {voucher.first().code} has been updated")
        return redirect(reverse("public_interface:browse"))

    if voucher:
        context["update"] = True
    else:
        context["update"] = False
    context["form"] = form
    return render(request, 'admin_interface/add_voucher.html', context)


def clean_post_data(request: HttpRequest) -> Dict[str, Any]:
    """Prepares request.POST data to be used to update a Vouchers object"""
    data = dict()
    for key, value in request.POST.items():
        if key in ["latitude", "longitude"]:
            try:
                value = float(value)
            except ValueError:
                value = None
        if key in ["min_altitude", "max_altitude"]:
            try:
                value = float(value)
            except ValueError:
                value = None
        if key in ["date_extraction"]:
            try:
                value = datetime.strptime(value, "%Y-%m-%d")
            except ValueError:
                value = None
        if key not in ["csrfmiddlewaretoken", "update"]:
            data[key] = value
    return data
