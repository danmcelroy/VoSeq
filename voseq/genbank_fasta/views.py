import json

from django.shortcuts import render
from django.http import HttpResponseRedirect

from core.utils import get_version_stats
from .forms import GenBankFastaForm
from . import utils


def index(request):
    version, stats = get_version_stats()

    form = GenBankFastaForm()
    return render(request,
                  'genbank_fasta/index.html',
                  {
                      'form': form,
                      'version': version,
                      'stats': stats,
                  },
                  )


def results(request):
    version, stats = get_version_stats()

    if request.method == 'POST':
        form = GenBankFastaForm(request.POST)

        if form.is_valid():
            cleaned_data = form.cleaned_data
            print(cleaned_data)

            voucher_codes = utils.get_voucher_codes(cleaned_data)
            gene_codes = utils.get_gene_codes(cleaned_data)

            print(voucher_codes, gene_codes)
            res = utils.Results(voucher_codes, gene_codes)
            res.get_datasets()
            items_with_accession = res.items_with_accession
            fasta = res.fasta
            protein = res.protein

            return render(request, 'genbank_fasta/results.html',
                          {
                              'items_with_accession': items_with_accession,
                              'fasta': fasta,
                              'protein': protein,
                              'version': version,
                              'stats': stats,
                          },
                          )
        else:
            return render(request, 'genbank_fasta/index.html',
                          {
                              'form': form,
                              'version': version,
                              'stats': stats,
                          },
                          )

    return HttpResponseRedirect('/genbank_fasta/')
