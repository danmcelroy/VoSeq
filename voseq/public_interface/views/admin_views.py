from datetime import datetime
import logging
from typing import Dict, Any

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpRequest
from django.shortcuts import render, redirect, reverse

from core.utils import get_context
from public_interface.forms.admin_forms import VoucherForm, SequenceForm
from public_interface.models import Vouchers, Sequences


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


@login_required
def add_sequence(request: HttpRequest, sequence_id: str = None) -> HttpResponse:
    context = get_context(request)
    form = None
    sequence = None
    if sequence_id:
        sequence = Sequences.objects.filter(id=sequence_id, user=request.user)
        if sequence.exists():
            form = SequenceForm(
                request.POST or None,
                initial=sequence.values()[0],
                user=request.user,
            )

    if not form:
        form = SequenceForm(request.POST or None, user=request.user)

    if form.is_valid():
        sequence = form.save()
        sequence.user = request.user
        sequence.save()
        messages.add_message(request, messages.INFO,
                             f"The sequence of code {sequence.code} "
                             f"{sequence.gene_code} has been created")
        return redirect(reverse("public_interface:browse"))
    elif "update" in request.POST:
        # we are editing this voucher
        data = clean_post_data(request)
        sequence = Sequences.objects.filter(
            code__code=request.POST.get("code"),
            gene_code=request.POST.get("gene_code"),
        )
        sequence.update(**data)
        messages.add_message(request, messages.INFO,
                             f"The sequence of code {sequence.first().code} "
                             f"{sequence.first().gene_code} has been updated")
        return redirect(reverse("public_interface:browse"))
    else:
        print(form.errors)

    if sequence:
        context["update"] = True
    else:
        context["update"] = False
    context["form"] = form
    context["sequence"] = sequence
    return render(request, 'admin_interface/add_sequence.html', context)


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
