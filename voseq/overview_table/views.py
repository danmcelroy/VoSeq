from django.shortcuts import render

from core.utils import get_version_stats
from core.utils import get_username
from public_interface.models import Vouchers


def index(request):
    version, stats = get_version_stats()
    username = get_username(request)

    return render(request,
                  'overview_table/index.html',
                  {
                      'username': username,
                      'version': version,
                      'stats': stats,
                  },
                  )
