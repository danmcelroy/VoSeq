import csv
import datetime
import json

from django.http import HttpResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt

from core.utils import get_version_stats
from public_interface.models import Vouchers


def index(request):
    version, stats = get_version_stats()
    return render(request, 'gbif/index.html',
                  {
                      'version': version,
                      'stats': stats,
                  },
                  )


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


def get_data_count():
    voucher_count = Vouchers.objects.count()
    return voucher_count


def create_excel_file():
    today = datetime.date.today()
    filename = 'data_for_GBIF_' + datetime.datetime.strftime(today, '%Y%m%d') + '.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="' + filename + '"'

    writer = csv.writer(response)

    row = ['Code', 'Order', 'Superfamily', 'Family', 'Subfamily', 'Tribe', 'Subtribe', 'Genus',
           'Species', 'Subspecies', 'TypeSpecies', 'Country', 'Specific_Locality', 'Latitude',
           'Longitude', 'Max. Altitude', 'Min. Altitude', 'Collector', 'Date_of_Collection',
           'Voucher_Locality', 'Host_Organism', 'Sex', 'Voucher_State', 'Voucher_code_from_others',
           'Code from BOLD',
           'Date_of_DNA_extraction', 'Extractor', 'Extraction #', 'Extraction_Vial',
           'Published_in', 'Notes']

    writer.writerow(row)

    vouchers = Vouchers.objects.all()
    for voucher in vouchers:
        row = [voucher.code]
        row.append(voucher.orden)
        row.append(voucher.superfamily)
        row.append(voucher.family)
        row.append(voucher.subfamily)
        row.append(voucher.tribe)
        row.append(voucher.subtribe)
        row.append(voucher.genus)
        row.append(voucher.species)
        row.append(voucher.subspecies)
        row.append(voucher.typeSpecies)
        row.append(voucher.country)
        row.append(voucher.specificLocality)
        row.append(voucher.latitude)
        row.append(voucher.longitude)
        row.append(voucher.max_altitude)
        row.append(voucher.min_altitude)
        row.append(voucher.collector)
        row.append(voucher.dateCollection)
        row.append(voucher.voucherLocality)
        row.append(voucher.hostorg)
        row.append(voucher.sex)
        row.append(voucher.voucher)
        row.append(voucher.voucherCode)
        row.append(voucher.code_bold)
        row.append(voucher.dateExtraction)
        row.append(voucher.extractor)
        row.append(voucher.extractionTube)
        row.append(voucher.publishedIn)
        row.append(voucher.notes)
        writer.writerow(row)
    return response


def results(request):
    return ''
