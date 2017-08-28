"""
Django settings for voseq project.

For more information on this file, see
https://docs.djangoproject.com/en/1.7/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.7/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os

BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DEBUG = True

MANAGERS = ()
ADMINS = ()

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.7/howto/deployment/checklist/


ALLOWED_HOSTS = []

# for testing in Travis CI
TRAVIS = False

# Application definition

INSTALLED_APPS = (
    # added
    'suit',
    'haystack',
    'crispy_forms',

    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',
    'django.contrib.sites',

    # my apps
    'core',
    'public_interface',
    'create_dataset',
    'blast_local',
    'blast_local_full',
    'blast_ncbi',
    'blast_new',
    'stats',
    'view_genes',
    'genbank_fasta',
    'gene_table',
    'voucher_table',
    'gbif',
    'overview_table',

    'registration',
    'easy_thumbnails',
    'import_export',

)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'debug_toolbar.middleware.DebugToolbarMiddleware',
)

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            'core/templates',
            'voseq/public_interface/templates/public_interface',
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
            'debug': False,
        }
    },
]
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'file': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'filename': '/tmp/voseq.log',
        },
    },
    'loggers': {
        'django': {
            'handlers': ['file'],
            'level': 'DEBUG',
            'propagate': True,
        },
    },
}
ROOT_URLCONF = 'voseq.urls'

WSGI_APPLICATION = 'voseq.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.7/ref/settings/#databases

# Internationalization
# https://docs.djangoproject.com/en/1.7/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

CRISPY_TEMPLATE_PACK = 'bootstrap3'


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.7/howto/static-files/
MEDIA_URL = '/media/'

STATIC_URL = '/static/'
STATICFILES_DIRS = (
    BASE_DIR,
    os.path.join(BASE_DIR, '..', 'public_interface'),
    os.path.join(BASE_DIR, '..', 'create_dataset'),
)

THUMBNAIL_ALIASES = {
    '': {
        'thumb': {'size': (200, 200), 'crop': False},
    },
}
# Get your API key from here:
# https://developers.google.com/maps/documentation/javascript/tutorial#api_key
# so that you can show Google Maps in your voucher pages.
GOOGLE_MAPS_API_KEY = 'fake api key'

FLICKR_API_KEY = "fake api key"
FLICKR_API_SECRET = "fake api secret"


# This VoSeq version
def get_version():
    import voseq
    return voseq.__version__

VERSION = get_version()


TESTING = False

# Django registration redux
ACCOUNT_ACTIVATION_DAYS = 7  # One-week activation window; you may, of course, use a different value.
REGISTRATION_AUTO_LOGIN = True  # Automatically log the user in.
SITE_ID = 1
LOGIN_REDIRECT_URL = '/browse/'

# Change this after obtaining VoSeq and before deployments to a production server
SECRET_KEY = 'test_key'
