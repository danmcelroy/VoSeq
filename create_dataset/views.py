import logging
import os
import re
from uuid import uuid4

from celery import chord
from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.http import HttpResponse

from core.utils import get_context
from public_interface.tasks import log_email_error, notify_user
from .forms import CreateDatasetForm
from .models import Dataset
from .tasks import create_dataset


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
            task_id = schedule_dataset(form.cleaned_data, request.user)
            context['task_id'] = task_id
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


def schedule_dataset(cleaned_data, user):
    taxonset_id = cleaned_data['taxonset'].id
    geneset_id = cleaned_data['geneset'].id
    gene_codes_ids = list(cleaned_data['gene_codes'].values_list('id', flat=True))
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

    dataset_obj = Dataset.objects.create(
        user=user
    )

    tasks = chord(
        header=create_dataset.si(
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
            dataset_obj.id,
        ).on_error(log_email_error.s(user.id)),
        body=notify_user.si(task_id, user.id)
    )
    tasks.apply_async(task_id=task_id)
    return task_id
