import os
import re

from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt

from core.utils import get_version_stats
from core.utils import get_gene_codes
from core.utils import get_voucher_codes
from .forms import GenBankFastaForm
from create_dataset.utils import CreateDataset
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


@csrf_exempt
def results(request):
    version, stats = get_version_stats()

    if request.method == 'POST':
        form = GenBankFastaForm(request.POST)

        if form.is_valid():
            cleaned_data = form.cleaned_data
            cleaned_data['file_format'] = 'GenbankFASTA'
            cleaned_data['number_genes'] = ''
            cleaned_data['aminoacids'] = True
            cleaned_data['positions'] = 'ALL'
            cleaned_data['partition_by_positions'] = 'ONE'
            cleaned_data['taxon_names'] = ['CODE', 'GENUS', 'SPECIES']
            cleaned_data['outgroup'] = ''

            print(">>>>", cleaned_data)
            dataset_creator = CreateDataset(cleaned_data)
            dataset = dataset_creator.dataset_str
            errors = dataset_creator.errors
            warnings = dataset_creator.warnings

            dataset_file_abs = dataset_creator.dataset_file
            if dataset_file_abs is not None:
                dataset_file = re.search('([A-Z]+_[a-z0-9]+\.txt)', dataset_file_abs).groups()[0]
            else:
                dataset_file = False

            return render(request, 'genbank_fasta/results.html',
                          {
                              'items_with_accession': '',
                              'fasta': '',
                              'fasta_file': dataset_file,
                              'protein': '',
                              'protein_file': '',
                              'warnings': '',
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


def serve_file(request, file_name):
    cwd = os.path.dirname(__file__)
    fasta_file = os.path.join(cwd,
                              'fasta_files',
                              file_name,
                              )
    response = HttpResponse(open(fasta_file, 'r').read(), content_type='application/text')
    response['Content-Disposition'] = 'attachment; filename=voseq_genbank.fasta'

    if os.path.isfile(fasta_file):
        os.remove(fasta_file)
    return response
