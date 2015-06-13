#!/bin/bash

# Script to set up dependencies for Django on Vagrant.

PGSQL_VERSION=9.3

# Need to fix locale so that Postgres creates databases in UTF-8
cp -p /vagrant_data/etc-bash.bashrc /etc/bash.bashrc
locale-gen en_GB.UTF-8
dpkg-reconfigure locales

export LANGUAGE=en_GB.UTF-8
export LANG=en_GB.UTF-8
export LC_ALL=en_GB.UTF-8

# Install essential packages from Apt
apt-get update -y
# Python dev packages
apt-get install -y build-essential python python3-dev python-setuptools python-pip
# Dependencies for image processing with Pillow (drop-in replacement for PIL)
# supporting: jpeg, tiff, png, freetype, littlecms
apt-get install -y libjpeg-dev libtiff-dev zlib1g-dev libfreetype6-dev liblcms2-dev
# Git (we'd rather avoid people keeping credentials for git commits in the repo, but sometimes we need it for pip requirements that aren't in PyPI)
apt-get install -y git

# Postgresql
if ! command -v psql; then
    apt-get install -y postgresql-$PGSQL_VERSION libpq-dev
    cp /vagrant_data/pg_hba.conf /etc/postgresql/$PGSQL_VERSION/main/
    /etc/init.d/postgresql reload
fi

# virtualenv global setup
if ! command -v pip; then
    easy_install -U pip
fi
if [[ ! -f /usr/local/bin/virtualenv ]]; then
    easy_install virtualenv virtualenvwrapper stevedore virtualenv-clone
fi

# bash environment global setup
cp -p /vagrant_data/bashrc /home/vagrant/.bashrc

# Cleanup
apt-get clean


# Virtualenv for VoSeq
if [[ ! -e /home/vagrant/.virtualenvs/voseq ]]; then
    su - vagrant -c "source /usr/local/bin/virtualenvwrapper.sh &&            \
        mkvirtualenv -p /usr/bin/python3 voseq && pip install pip --upgrade && \
        pip install -r /vagrant/requirements/base.txt"
fi

apt-get -y install openjdk-7-jdk openjdk-7-jre
if [[ ! -f /etc/init.d/elasticsearch ]]; then
    wget https://download.elastic.co/elasticsearch/elasticsearch/elasticsearch-1.6.0.deb && \
        dpkg -i elasticsearch-1.6.0.deb && /etc/init.d/elasticsearch start
fi

