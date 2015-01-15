from django.shortcuts import render
from django.conf import settings

from public_interface.models import Genes
from stats.models import Stats


def index(request):
    VERSION = settings.VERSION
    try:
        STATS = Stats.objects.get(pk=1)
    except Stats.DoesNotExist:
        STATS = ''

    queryset = Genes.objects.all()
    print(queryset.values())
    return render(request, 'view_genes/index.html',
                  {
                      'result': queryset,
                      'version': VERSION,
                      'stats': STATS,
                  },
                  )
