import json

from django.http import HttpResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt

from core.utils import get_context
from .utils import get_data_count
from .utils import create_excel_file


def index(request):
    context = get_context(request)
    return render(request, 'gbif/index.html', context)


@csrf_exempt
def dump_data(request):
    try:
        wanted = request.GET['request']
    except KeyError:
        msg = json.dumps({'result': 'error'})
        return HttpResponse(msg, content_type='application/json')

    if wanted == 'count_data':
        the_data = get_data_count()
        msg = json.dumps({
            'result': True,
            'count': the_data,
        })
        return HttpResponse(msg, content_type='application/json')
    if wanted == 'make_file':
        response = create_excel_file()
        return response
