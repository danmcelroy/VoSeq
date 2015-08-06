# -*- coding: utf-8 -*-
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import voseq

version = voseq.__version__

requirements = [
    'biopython',
    'Django',
    'pyprind',
    'elasticsearch',
    'Unipath',
    'psycopg2',
    'dataset',
    'django-haystack',
    'django-debug-toolbar',
    'pytz',
    'django-suit',
]

test_requirements = [
    'coverage',
    'nose',
    'pep8',
    'coveralls',
    'Sphinx',
]

setup(
    name='VoSeq',
    version=version,
    author='Carlos Pena, Tobias Malm, Victor Mezarino',
    author_email='carlosp420@gmail.com',
    packages=[
        'voseq',
    ],
    include_package_data=True,
    install_requires=requirements,
    license='BSD',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    zip_safe=False,
    scripts=['voseq/manage.py'],
    test_require=test_requirements,
)
