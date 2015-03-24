import re

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
        form = CreateDatasetForm(request.POST)

        if form.is_valid():
            print(">>>>", form.cleaned_data)
            dataset_creator = CreateDataset(form.cleaned_data)
            dataset = dataset_creator.dataset_str
            errors = dataset_creator.errors
            warnings = dataset_creator.warnings

            phylip_file = dataset_creator.phylip_partition_file
            if phylip_file is not None:
                phylip_partition_file = re.search('(phylip_[a-z0-9]+_partitions\.phy)', phylip_file).groups()[0]

            return render(request, 'create_dataset/results.html',
                          {
                              'phylip_partitions_file': phylip_partition_file,
                              'dataset': dataset,
                              'errors': errors,
                              'warnings': warnings,
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
