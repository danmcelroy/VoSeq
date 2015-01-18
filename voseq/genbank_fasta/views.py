from django.shortcuts import render
from django.http import HttpResponseRedirect

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


def results(request):
    version, stats = get_version_stats()

    if request.method == 'POST':
        form = GenBankFastaForm(request.POST)

        if form.is_valid():
            cleaned_data = form.cleaned_data
            print(cleaned_data)

            result = ''

            return render(request, 'genbank_fasta/results.html',
                          {
                              'result': result,
                              'version': version,
                              'stats': stats,
                          },
                          )
        else:
            print(form)
            return render(request, 'genbank_fasta/index.html',
                          {
                              'form': form,
                              'version': version,
                              'stats': stats,
                          },
                          )

    return HttpResponseRedirect('/genbank_fasta/')
