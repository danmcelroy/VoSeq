import logging
import os
import re
from uuid import uuid4

from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.http import HttpResponse

from core.utils import get_context
from .forms import CreateDatasetForm
from .tasks import create_dataset
from .utils import CreateDataset


log = logging.getLogger(__name__)


@login_required
def index(request):
    form = CreateDatasetForm()
    context = get_context(request)
    context["form"] = form

    return render(request, 'create_dataset/index.html', context)


@login_required
def results(request):
    context = get_context(request)

    if request.method == 'POST':
        form = CreateDatasetForm(request.POST)

        if form.is_valid():
            dataset_format = form.cleaned_data['file_format']
            job_id = schedule_dataset(form.cleaned_data)
            dataset_creator = CreateDataset(form.cleaned_data)
            dataset = "{}{}{}".format(
                dataset_creator.dataset_str[0:1500],
                '\n...\n\n\n',
                '#######\nComplete dataset file available for download.\n#######',
            )
            errors = dataset_creator.errors
            warnings = set(dataset_creator.warnings)

            dataset_file_abs = dataset_creator.dataset_file
            if dataset_file_abs is not None:
                dataset_file = re.search(
                    '([A-Z]+_[a-z0-9]+\.txt)',
                    dataset_file_abs
                ).groups()[0]
            else:
                dataset_file = False

            context['dataset_file'] = dataset_file
            context['charset_block'] = dataset_creator.charset_block
            context['dataset'] = dataset
            context['dataset_format'] = dataset_format
            context['errors'] = errors
            context['warnings'] = warnings
            return render(request, 'create_dataset/results.html', context)
        else:
            log.debug("invalid form")
            context["form"] = form
            return render(request, 'create_dataset/index.html', context)
    else:
        return HttpResponseRedirect('/create_dataset/')


@login_required
def serve_file(request, file_name):
    final_name = guess_file_extension(file_name)
    cwd = os.path.dirname(__file__)
    dataset_file = os.path.join(cwd,
                                'dataset_files',
                                file_name,
                                )
    if os.path.isfile(dataset_file):
        response = HttpResponse(open(dataset_file, 'r').read(), content_type='application/text')
        response['Content-Disposition'] = 'attachment; filename={}'.format(final_name)
        os.remove(dataset_file)
        return response
    else:
        return render(request, 'create_dataset/missing_file.html')


def guess_file_extension(file_name):
    try:
        prefix = re.search('^(\w+)\_', file_name).group()
    except AttributeError:
        return file_name

    if prefix == 'MEGA_':
        extension = 'meg'
    else:
        return file_name

    name = file_name.replace('.txt', '')
    return '{0}.{1}'.format(name, extension)


def schedule_dataset(cleaned_data):
    taxonset_id = cleaned_data['taxonset'].id
    geneset_id = cleaned_data['geneset'].id
    gene_codes_ids = cleaned_data['gene_codes'].values_list('id', flat=True)
    voucher_codes = cleaned_data['voucher_codes']
    file_format = cleaned_data['file_format']
    outgroup = cleaned_data['outgroup']
    positions = cleaned_data['positions']
    partition_by_positions = cleaned_data['partition_by_positions']
    translations = cleaned_data['translations']
    aminoacids = cleaned_data['aminoacids']
    degen_translations = cleaned_data['degen_translations']
    special = cleaned_data['special']
    taxon_names = cleaned_data['taxon_names']
    number_genes = cleaned_data['number_genes']
    introns = cleaned_data['introns']
    task_id = str(uuid4())

    task = create_dataset.si(
        (
            taxonset_id,
            geneset_id,
            gene_codes_ids,
            voucher_codes,
            file_format,
            outgroup,
            positions,
            partition_by_positions,
            translations,
            aminoacids,
            degen_translations,
            special,
            taxon_names,
            number_genes,
            introns,
        ),
        task_id=task_id,
    )
    task.apply_async()
    return task_id
