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


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.7/howto/deployment/checklist/



###################### FOR LOCAL DB CUSTOMIZATION #####################

import json
from django.core.exceptions import ImproperlyConfigured

SECRET_FILE = os.path.join(BASE_DIR, '..', 'config.json')

with open(SECRET_FILE) as f:
    secrets = json.loads(f.read())


def get_secret(setting, secrets=secrets):
    try:
        return secrets[setting]
    except KeyError:
        error_msg = "Set the {0} environment variable".format(setting)
        raise ImproperlyConfigured(error_msg)

#######################################################################


# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = get_secret("SECRET_KEY")

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

TEMPLATE_DIRS = (
    os.path.join(BASE_DIR, 'templates'),
)

ALLOWED_HOSTS = []


# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',

    # added
    'haystack',
    'debug_toolbar',

    # my apps
    'public_interface',
    'core',
    'stats',
    'blast_local',
    'blast_local_full',
    'blast_ncbi',
    'blast_new',
    'view_genes',
    'genbank_fasta',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.contrib.messages.context_processors.messages',
    'django.contrib.auth.context_processors.auth',
)

ROOT_URLCONF = 'voseq.urls'

WSGI_APPLICATION = 'voseq.wsgi.application'

HAYSTACK_CONNECTIONS = {
    'default': {
        'ENGINE': 'haystack.backends.elasticsearch_backend.ElasticsearchSearchEngine',
        'URL': 'http://127.0.0.1:9200/',
        'INDEX_NAME': 'haystack',
    }
}

HAYSTACK_DEFAULT_OPERATOR = 'AND'


# Database
# https://docs.djangoproject.com/en/1.7/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': get_secret('DB_NAME'),
        'USER': get_secret('DB_USER'),
        'PASSWORD': get_secret('DB_PASS'),
        'HOST': get_secret('DB_HOST'),
        'PORT': get_secret('DB_PORT'),
    }
}

# Internationalization
# https://docs.djangoproject.com/en/1.7/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.7/howto/static-files/

STATIC_URL = '/static/'
STATICFILES_DIRS = (
    BASE_DIR,
)

# Get your Google API key so that you can show Google Maps in your voucher pages
# https://developers.google.com/maps/documentation/javascript/tutorial#api_key

GOOGLE_MAPS_API_KEY = get_secret("GOOGLE_MAPS_API_KEY")


# This VoSeq version
VERSION = '2.0.0'




