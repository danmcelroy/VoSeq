from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponseRedirect

from .utils import BLASTNew
from stats.models import Stats
from .forms import BLASTNewForm


def index(request):
    VERSION = settings.VERSION
    try:
        STATS = Stats.objects.get(pk=1)
    except Stats.DoesNotExist:
        STATS = ''

    form = BLASTNewForm()
    return render(request, 'blast_new/index.html',
                  {
                      'form': form,
                      'version': VERSION,
                      'stats': STATS,
                  },
                  )


def results(request):
    VERSION = settings.VERSION
    try:
        STATS = Stats.objects.get(pk=1)
    except Stats.DoesNotExist:
        STATS = ''

    if request.method == 'POST':
        form = BLASTNewForm(request.POST)

        if form.is_valid():
            cleaned_data = form.cleaned_data
            print(cleaned_data)

            blast = BLASTNew('new', cleaned_data['name'], cleaned_data['sequence'],
                             cleaned_data['gene_codes'])
            blast.save_seqs_to_file()

            if blast.is_blast_db_up_to_date() is False:
                blast.create_blast_db()

            blast.save_query_to_file()
            blast.do_blast()
            result = blast.parse_blast_output()
            blast.delete_query_output_files()
            return render(request, 'blast_new/results.html',
                          {
                              'result': result,
                              'version': VERSION,
                              'stats': STATS,
                          },
                          )
        else:
            return render(request, 'blast_new/index.html',
                          {
                              'version': VERSION,
                              'stats': STATS,
                              'form': form,
                          },
                          )

    return HttpResponseRedirect('/blast_new/')
