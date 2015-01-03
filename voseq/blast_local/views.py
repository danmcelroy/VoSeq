

from django.shortcuts import render


def index(request, voucher_code, gene_code):
    return render(request, 'public_interface/base.html')
