from django.shortcuts import render

from core.utils import get_version_stats


def index(request):
    version, stats = get_version_stats()
    return render(request, 'gbif/index.html',
                  {
                      'version': version,
                      'stats': stats,
                  },
                  )


def results(request):
    return ''
