from django.db.models import Q


class IssueSearch(object):
    """
    Class to do an advanced search.
    Taken from http://www.slideshare.net/tpherndon/django-search-presentation
    """
    def __init__(self, search_data):
        self.__dict__.update(search_data)

    def search_keywrods(self, q):
        if self.keywords:
            words = self.keywords.split()
            code_q = Q()
            specific_locality_q = Q()

            for word in words:
                code_q = code_q | Q(title__icontains=word)
                specific_locality_q = specific_locality_q | Q(specific_locality__icontains=word)
            keyword_q = code_q | specific_locality_q
            q = q & keyword_q
        return q
