"""
To be used in Travis CI so we can skip running the tests for NCBI blast, which
fails most of the time due to network problems.
"""

import sys

from .base import *


print('Testing in Travis')
TESTING = len(sys.argv) > 1 and sys.argv[1] == 'test'
TRAVIS = True

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': 'test.db',
    }
}
