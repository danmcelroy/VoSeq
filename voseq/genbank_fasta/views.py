import os
import re

from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt

from core.utils import get_version_stats
from core.utils import get_username
from .forms import GenBankFastaForm
from create_dataset.utils import CreateDataset


@login_required
def index(request):
    version, stats = get_version_stats()
    username = get_username(request)

    form = GenBankFastaForm()
    return render(request,
                  'genbank_fasta/index.html',
                  {
                      'username': username,
                      'form': form,
                      'version': version,
                      'stats': stats,
                  },
                  )


@login_required
@csrf_exempt
def results(request):
    version, stats = get_version_stats()
    username = get_username(request)

    if request.method == 'POST':
        form = GenBankFastaForm(request.POST)

        if form.is_valid():
            cleaned_data = form.cleaned_data
            cleaned_data['file_format'] = 'GenBankFASTA'
            cleaned_data['number_genes'] = ''
            cleaned_data['translations'] = False
            cleaned_data['aminoacids'] = False
            cleaned_data['positions'] = 'ALL'
            cleaned_data['partition_by_positions'] = 'by gene'
            cleaned_data['taxon_names'] = ['CODE', 'GENUS', 'SPECIES']
            cleaned_data['outgroup'] = ''

            dataset_creator = CreateDataset(cleaned_data)
            dataset = dataset_creator.dataset_str
            dataset_short = dataset[0:1500] + '\n...\n\n\n' + '#######\nComplete dataset file available for download.\n#######'
            errors = dataset_creator.errors
            warnings = dataset_creator.warnings
            dataset_file_abs = dataset_creator.dataset_file
            if dataset_file_abs is not None:
                dataset_file = os.path.basename(dataset_file_abs)
            else:
                dataset_file = False

            cleaned_data['aminoacids'] = True
            dataset_creator = CreateDataset(cleaned_data)
            aa_dataset = dataset_creator.dataset_str
            aa_dataset_file_abs = dataset_creator.dataset_file
            if aa_dataset_file_abs is not None:
                aa_dataset_file = os.path.basename(aa_dataset_file_abs)
            else:
                aa_dataset_file = False

            return render(request, 'genbank_fasta/results.html',
                          {
                              'username': username,
                              'items_with_accession': '',
                              'dataset': dataset_short,
                              'fasta_file': dataset_file,
                              'protein': aa_dataset,
                              'errors': errors,
                              'protein_file': aa_dataset_file,
                              'warnings': warnings,
                              'version': version,
                              'stats': stats,
                          },
                          )
        else:
            return render(request, 'genbank_fasta/index.html',
                          {
                              'username': username,
                              'form': form,
                              'version': version,
                              'stats': stats,
                          },
                          )

    return HttpResponseRedirect('/genbank_fasta/')


@login_required
def serve_file(request, file_name):
    cwd = os.path.dirname(__file__)
    fasta_file = os.path.join(cwd,
                              '..',
                              'create_dataset',
                              'dataset_files',
                              file_name,
                              )
    response = HttpResponse(open(fasta_file, 'r').read(), content_type='application/text')
    response['Content-Disposition'] = 'attachment; filename=voseq_genbank.fasta'

    if os.path.isfile(fasta_file):
        os.remove(fasta_file)
    return response
