from itertools import chain
import json
import logging

from django.contrib.auth.decorators import login_required
from django.db.models import Q
from django.http import HttpResponseRedirect, Http404, HttpResponse
from django.views.decorators.csrf import csrf_protect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.shortcuts import render, redirect
from django.conf import settings

from haystack.forms import SearchForm
from haystack.query import ValuesSearchQuerySet

from core.utils import get_context
from public_interface.utils import get_simple_query, get_correct_url_query, get_voucher_code_list
from public_interface.models import Vouchers, FlickrImages, LocalImages, Sequences, Primers
from public_interface.forms.admin_forms import VoucherForm


log = logging.getLogger(__name__)


@login_required
def add_voucher(request):
    context = get_context(request)
    form = VoucherForm(request.POST or None)

    context["form"] = form
    return render(request, 'admin_interface/add_voucher.html', context)
