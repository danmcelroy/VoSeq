from datetime import datetime
import logging
from typing import Dict, Any

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpRequest
from django.shortcuts import render, redirect, reverse

from core.utils import get_context
from public_interface.forms.admin_forms import VoucherForm, SequenceForm, GeneForm
from public_interface.models import Vouchers, Sequences, Genes


log = logging.getLogger(__name__)


@login_required
def add_gene(request: HttpRequest, gene_id: str = None) -> HttpResponse:
    context = get_context(request)
    form = None
    gene = None
    if gene_id:
        gene = Genes.objects.filter(id=gene_id)
        if gene:
            form = GeneForm(request.POST or None, initial=gene.values()[0])

    if not form:
        form = GeneForm(request.POST or None)

    if form.is_valid():
        voucher = form.save()
        voucher.user = request.user
        voucher.save()
        return HttpResponse("form is valid")
    elif "update" in request.POST:
        data = clean_post_data(request)
        gene = Genes.objects.filter(code=request.POST.get("code"))
        gene.update(**data)
        messages.add_message(request, messages.INFO,
                             f"The gene of code {gene.first().gene_code} has been updated")
        return redirect(reverse("public_interface:browse"))

    if gene:
        context["update"] = True
    else:
        context["update"] = False
    context["form"] = form
    return render(request, 'admin_interface/add_gene.html', context)

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

    if "update" in request.POST:
        # we are editing this voucher
        data = clean_post_data(request)
        sequence = Sequences.objects.filter(id=data.pop("sequence_id"))
        sequence.update(**data)
        messages.add_message(request, messages.INFO,
                             f"The sequence of code {sequence.first().code} "
                             f"{sequence.first().gene_code} has been updated")
        return redirect(reverse("public_interface:browse"))
    elif form.is_valid():
        sequence = form.save()
        sequence.user = request.user
        sequence.save()
        messages.add_message(request, messages.INFO,
                             f"The sequence of code {sequence.code} "
                             f"{sequence.gene_code} has been created")
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
        if key == "gene_code":
            value = Genes.objects.get(id=value).gene_code

        if value == "None":
            value = None
        elif value == "True":
            value = True
        elif value == "False":
            value = False

        if key not in ["csrfmiddlewaretoken", "update"]:
            data[key] = value
    return data
