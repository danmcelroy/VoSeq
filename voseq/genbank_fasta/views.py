from django.shortcuts import render

from core.utils import get_version_stats
from .forms import GenBankFastaForm


def index(request):
    version, stats = get_version_stats()

    form = GenBankFastaForm()
    return render(request,
                  'genbank_fasta/index.html',
                  {
                      'form': form,
                      'version': version,
                      'stats': stats,
                  },
                  )
