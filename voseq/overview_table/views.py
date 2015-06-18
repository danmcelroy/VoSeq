from django.shortcuts import render

from core.utils import get_version_stats
from core.utils import get_username
from public_interface.models import Genes
from overview_table.models import OverviewTable


def index(request):
    version, stats = get_version_stats()
    username = get_username(request)

    genes = Genes.objects.all().order_by('gene_code')
    vouchers = OverviewTable.objects.all()
    return render(request,
                  'overview_table/index.html',
                  {
                      'username': username,
                      'version': version,
                      'stats': stats,
                      'data': vouchers,
                      'genes': genes,
                  },
                  )
