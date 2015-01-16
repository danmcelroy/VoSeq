from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponseRedirect

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
    if request.method == 'POST':
        form = BLASTNewForm(request.POST)
        if form.is_valid():
            print(request.POST)
            return render(request,
                          'blast_new/results.html')

    return HttpResponseRedirect('/blast_new/')
