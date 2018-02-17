import logging

from django.shortcuts import render
from django.http import HttpResponseRedirect

from core.utils import get_version_stats, get_username
from .utils import BLASTNew
from .forms import BLASTNewForm


log = logging.getLogger(__name__)


def index(request):
    version, stats = get_version_stats()
    username = get_username(request)

    form = BLASTNewForm()
    return render(request, 'blast_new/index.html',
                  {
                      'username': username,
                      'form': form,
                      'version': version,
                      'stats': stats,
                  },
                  )


def results(request):
    version, stats = get_version_stats()
    username = get_username(request)

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
            return render(request, 'blast_new/results.html',
                          {
                              'username': username,
                              'result': result,
                              'version': version,
                              'stats': stats,
                          },
                          )
        else:
            return render(request, 'blast_new/index.html',
                          {
                              'username': username,
                              'form': form,
                              'version': version,
                              'stats': stats,
                          },
                          )

    return HttpResponseRedirect('/blast_new/')
