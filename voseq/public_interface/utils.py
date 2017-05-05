import re
from urllib.parse import parse_qs

from haystack.views import SearchView

from core.utils import get_version_stats
from core.utils import get_username


def get_simple_query(request):
    url_encoded_query = get_correct_url_query(request.GET.urlencode())
    simple_query = recover_keyword(url_encoded_query)
    return simple_query


def strip_page(url_encoded_query):
    this_query = re.sub('page=[0-2]+', '', url_encoded_query)
    this_query = this_query.replace('&&', '&')
    this_query = re.sub('^&', '', this_query)
    this_query = re.sub('&$', '', this_query)
    return this_query


def get_correct_url_query(url_encoded_query):
    this_query = strip_page(url_encoded_query)
    return this_query


def recover_keyword(url_encoded_query):
    simple_query = ""
    for key, value in parse_qs(url_encoded_query).items():
        if value[0] != "Select" and key != "page":
            simple_query += value[0] + " "
    return simple_query.strip()


def get_voucher_code_list(sqs):
    if sqs is None:
        return None
    try:
        # this is voucher queryset
        code = sqs[0].code
        code_list = "\n".join(sqs.values_list("code", flat=True))
    except TypeError:
        # this is sequences queryset
        code_list = "\n".join(sqs.values_list("code__code", flat=True))
    return code_list


class VoSeqSearchView(SearchView):
    def __init__(self, url_encoded_query, *args, **kwargs):
        self.searchqueryset = kwargs['searchqueryset']
        self.url_encoded_query = self.get_correct_url_query(url_encoded_query)
        self.simple_query = self.recover_keyword(url_encoded_query)
        self.voucher_code_list = self.get_voucher_code_list()
        super(VoSeqSearchView, self).__init__(*args, **kwargs)

    def get_correct_url_query(self, url_encoded_query):
        this_query = self.strip_page(url_encoded_query)
        return this_query

    def recover_keyword(self, url_encoded_query):
        this_query = re.sub('\w+=Select', ' ', url_encoded_query)
        this_query = re.sub('\w+=', ' ', this_query)
        this_query = this_query.replace('&', ' ')
        this_query = re.sub('\s+', ' ', this_query)
        return this_query.strip()

    def strip_page(self, url_encoded_query):
        this_query = re.sub('page=[0-2]+', '', url_encoded_query)
        this_query = this_query.replace('&&', '&')
        this_query = re.sub('^&', '', this_query)
        this_query = re.sub('&$', '', this_query)
        return this_query

    def get_voucher_code_list(self):
        if self.searchqueryset is None:
            return None
        code_list = ''
        for i in self.searchqueryset:
            code_list += i.code + '\n'
        return code_list

    def extra_context(self):
        version, stats = get_version_stats()
        username = get_username(self.request)
        return {
            'username': username,
            'voucher_code_list': self.voucher_code_list,
            'simple_query': self.simple_query,
            'url_encoded_query': self.url_encoded_query,
            'result_count': len(self.searchqueryset),
            'version': version,
            'stats': stats,
        }
