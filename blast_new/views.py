import logging

from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from django.http import HttpResponseRedirect

from core.utils import get_context
from .utils import BLASTNew
from .forms import BLASTNewForm


log = logging.getLogger(__name__)


@login_required
def index(request):
    form = BLASTNewForm()
    context = get_context(request)
    context["form"] = form
    return render(request, 'blast_new/index.html', context=context)


@login_required
def results(request):
    context = get_context(request)

    if request.method == 'POST':
        log.debug(request.POST)
        form = BLASTNewForm(request.POST)

        if form.is_valid():
            cleaned_data = form.cleaned_data

            blast = BLASTNew(
                blast_type='new',
                name=cleaned_data['name'],
                sequence=cleaned_data['sequence'],
                gene_codes=cleaned_data['gene_codes'],
            )
            blast.save_seqs_to_file()

            if not blast.is_blast_db_up_to_date():
                blast.create_blast_db()

            blast.save_query_to_file()
            blast.do_blast()
            result = blast.parse_blast_output()
            if not result:
                result = None
            blast.delete_query_output_files()
            context["result"] = result
            return render(request, 'blast_new/results.html', context)
        else:
            context["form"] = form
            return render(request, 'blast_new/index.html', context)

    return HttpResponseRedirect('/blast_new/')
