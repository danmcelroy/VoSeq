# -*- coding: utf-8 -*-
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import voseq

version = voseq.__version__

setup(
    name='VoSeq',
    version=version,
    author='Carlos Pena, Tobias Malm, Victor Mezarino',
    author_email='carlosp420@gmail.com',
    packages=[
        'voseq',
    ],
    include_package_data=True,
    install_requires=[
        'Django==1.7.3',
    ],
    zip_safe=False,
    scripts=['voseq/manage.py'],
)
