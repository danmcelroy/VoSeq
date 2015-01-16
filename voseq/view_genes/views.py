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
    return render(request, 'view_genes/index.html',
                  {
                      'result': queryset,
                      'version': VERSION,
                      'stats': STATS,
                  },
                  )


def gene(request, gene_code):
    VERSION = settings.VERSION
    try:
        STATS = Stats.objects.get(pk=1)
    except Stats.DoesNotExist:
        STATS = ''

    queryset = Genes.objects.filter(gene_code=gene_code)
    if len(queryset) < 1:
        item = ''
    else:
        item = queryset[0]

        if ';' in item.intron:
            introns = item.intron.split(';')
            j = 1
            out = ''
            for i in introns:
                out += '<b>Intron ' + str(j) + ':</b> ' + i + ' > '
                j += 1
            item.intron = out.rstrip(' > ')

    return render(request, 'view_genes/gene.html',
                  {
                      'item': item,
                      'version': VERSION,
                      'stats': STATS,
                  },
                  )
