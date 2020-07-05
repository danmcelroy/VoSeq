from django.http import JsonResponse

from public_interface.views import isearch_in_voucher


SEARCH_LIMIT = 20


def simple_search(request):
    query = request.GET.get('q')
    sqs = isearch_in_voucher(query)
    result = [str(voucher) for voucher in sqs[:SEARCH_LIMIT]]
    output = {'result': result}
    return JsonResponse(output)

