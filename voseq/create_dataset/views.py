from django.shortcuts import render
from django.http import HttpResponseRedirect

from core.utils import get_version_stats
from .forms import CreateDatasetForm
from .utils import CreateDataset


def index(request):
    form = CreateDatasetForm()

    return render(request,
                  'create_dataset/index.html',
                  {
                      'form': form,
                  },
                  )


def results(request):
    version, stats = get_version_stats()

    if request.method == 'POST':
        print(request.POST)
        form = CreateDatasetForm(request.POST)

        dataset_creator = CreateDataset()

        if form.is_valid():
            return render(request, 'create_dataset/results.html',
                          {
                              'dataset': dataset_creator.dataset_str,
                              'version': version,
                              'stats': stats,
                          },
                          )
        else:
            print("invalid form")
            return render(request, 'create_dataset/index.html',
                          {
                              'form': form,
                              'version': version,
                              'stats': stats,
                          },
                          )
    else:
        return HttpResponseRedirect('/create_dataset/')
