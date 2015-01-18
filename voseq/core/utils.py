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
