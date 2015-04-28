import sys

from .base import *


print('Testing')

SECRET_KEY = 'hola'

TESTING = len(sys.argv) > 1 and sys.argv[1] == 'test'

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': 'test.db',
    }
}
