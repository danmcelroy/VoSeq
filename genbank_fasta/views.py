import logging
import os

from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt

from core.utils import get_context
from .forms import GenBankFastaForm
from create_dataset.utils import CreateDataset


log = logging.getLogger(__name__)


@login_required
def index(request):
    form = GenBankFastaForm()
    context = get_context(request)
    context["form"] = form
    return render(request, 'genbank_fasta/index.html', context)


@login_required
@csrf_exempt
def results(request):
    context = get_context(request)

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
            dataset_short = dataset[0:1500] + '\n...\n\n\n' + '#######\nComplete dataset file available for download.\n#######'  # noqa
            errors = dataset_creator.errors
            warnings = dataset_creator.warnings
            dataset_file_abs = dataset_creator.dataset_file
            items_with_accession = dataset_creator.sequences_skipped
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

            context['items_with_accession'] = items_with_accession
            context['dataset'] = dataset_short
            context['fasta_file'] = dataset_file
            context['protein'] = aa_dataset
            context['errors'] = errors
            context['protein_file'] = aa_dataset_file
            context['warnings'] = warnings
            return render(request, 'genbank_fasta/results.html', context)
        else:
            context["form"] = form
            return render(request, 'genbank_fasta/index.html', context)

    return HttpResponseRedirect('/genbank_fasta/')


@login_required
def serve_file(request, file_name):
    log.debug("Requested file by user: {0}".format(request.user))
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
