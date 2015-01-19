import json

from django.shortcuts import render
from django.http import HttpResponseRedirect

from core.utils import get_version_stats
from .forms import GenBankFastaForm


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

            voucher_codes = []
            if cleaned_data['taxonset'] is not None:
                voucher_codes = json.loads(cleaned_data['taxonset'].taxonset_list)
            if cleaned_data['voucher_codes'] != '':
                voucher_codes += cleaned_data['voucher_codes'].splitlines()
            voucher_codes = set(voucher_codes)

            gene_codes = []
            if cleaned_data['geneset'] is not None:
                gene_codes = json.loads(cleaned_data['geneset'].geneset_list)
            if len(cleaned_data['gene_codes']) > 0:
                gene_codes += [i.gene_code for i in cleaned_data['gene_codes']]
            gene_codes = set(gene_codes)

            result = ''

            return render(request, 'genbank_fasta/results.html',
                          {
                              'result': result,
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
