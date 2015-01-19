import re

from django.conf import settings

from stats.models import Stats


def get_version_stats():
    """
    Returns version and database statistics for page footer.
    """
    version = settings.VERSION
    try:
        stats = Stats.objects.get(pk=1)
    except Stats.DoesNotExist:
        stats = ''

    return version, stats


def strip_question_marks(seq):
    """Having too many ambiguous characters will mess up DNA translation.
    """
    seq = re.sub('^\?+', '', seq)
    seq = re.sub('\?+$', '', seq)

    seq = re.sub('^N+', '', seq)
    seq = re.sub('N+$', '', seq)
    return seq
