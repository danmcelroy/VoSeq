|Waffle| |Dependency Status| |Coverage Status| |Landscape| |Docs|

+------------------+------------------+
| Windows          | Linux            |
+==================+==================+
| |Build status|   | |Build Status|   |
+------------------+------------------+

Contents
========

* `VoSeq is being rewritten`_
* `New features`_
* `Installation instructions`_
* `Test database for development`_
* `Start a test server`_
* `Migrate VoSeq database`_
* `Set-up a publicly available web server`_
* `Administrate the server`_


VoSeq is being rewritten
========================

We are rebuilding VoSeq from scratch. We decided to migrate from PHP to
Python by using the framework Django. We also moved from MySQL to
PostgreSQL.

You can still download the old VoSeq v1.7.4 from
`here <https://github.com/carlosp420/VoSeq/releases/tag/v1.7.4>`__. But
be aware that we will not be doing major maintenance of that code.

Here is a test installation of the old VoSeq (v1.7.0)
http://www.nymphalidae.net/VoSeq/

More details about the migration can be found in our `discussion
list <https://groups.google.com/forum/#!topic/voseq-discussion-list/wQ-E0Xcimgw>`__.

VoSeq 2.0.0 is the future!

New Features
============
Query suggestions for simple taxon searches:

.. image:: https://raw.githubusercontent.com/carlosp420/VoSeq/master/imgs/simple_search_suggestion.png

Installation instructions
=========================

These instructions assume that your libraries are up to date and that you have Python, pip, Java 7+ and virtual environments installed. Python3 is recommended.

**Step 1: get VoSeq.**
Clone or `download <https://github.com/carlosp420/VoSeq/releases>`__ VoSeq to your prefered directory.

**Step 2: create a virtual environment and install dependencies.**
To ensure that all the dependancies will work without conflict, it is best to install them within a virtual environment.

.. code:: shell

    mkvirtualenv -p /usr/bin/python voseq_environment
    cd /path/to/VoSeq
    workon voseq_environment
    pip install -r requirements/testing.txt

Exit the virtual environment for now to continue from the shell:

.. code:: shell

    deactivate

**Step 3: download and install elasticsearch.**
For elasticsearch, java needs to be installed. Mac users can download and install ``elasticsearch`` from here:
http://www.elasticsearch.org/overview/elkdownloads/ . In Linux, you can do:

.. code:: shell

    wget https://download.elastic.co/elasticsearch/elasticsearch/elasticsearch-1.5.2.deb
    sudo dpkg -i elasticsearch-1.5.2.deb

The bin directory of elasticsearch should be added automatically to your PATH. If not, add the following line to your ``.profile`` (Linux) or ``.bash_profile`` (macOSX) file:

.. code:: shell

    export PATH="$PATH:/parth/to/elasticsearch/bin/"

**Step 4: download, install and configure PostgreSQL.**
For macOSX users we recommend to do it by downloading the Postgres.app from http://postgresapp.com . Linux users can use apt-get:

.. code:: shell

    sudo apt-get install postgresql postgresql-contrib postgresql-server-dev-9.3

Create new role by typing:

.. code:: shell

    createuser --interactive

Enter the psql shell, create a password for this user and create a database for VoSeq:

.. code:: shell

    psql
    postgres=# ALTER ROLE postgres WITH PASSWORD 'hu8jmn3';
    postgres=# create database voseq;


In macOSX if you are using the Postgres.app, it my be enough to run:

.. code:: shell

    psql
    user.name=# CREATE DATABASE voseq;

To exit the psql shell:

.. code:: shell

    \q
    
Next, create a ``config.json`` file to keep the database variables:

.. code:: shell

    cd /path/to/Voseq
    touch config.json

and write in the following content:

.. code:: javascript

    {
    "SECRET_KEY": "create_a_secret_key",
    "DB_USER": "role_name",
    "DB_PASS": "create_a_database_password",
    "DB_NAME": "voseq",
    "DB_PORT": "5432",
    "DB_HOST": "localhost",
    "GOOGLE_MAPS_API_KEY": "get_a_google_map_api_key"
    }

If you followed the above instructions to the letter, the DB_USER will be "postgres" and the DB_PASS will be "hu8jmn3". It is of recommended to come up with your own password. Instructions to obtain a personal google map browser API key can be found `here <https://developers.google.com/maps/documentation/javascript/tutorial#api_key>`__. 

After following these four steps everything should be installed and ready to run. You can now choose to either continue with adding real data migrated from VoSeq 1.x and setting up a publicly available web server, or to first add some test data and test the set-up with a lightweight local server included in the VoSeq package.

Test database for development
=============================

You can use test data to populate your PostgreSQL database, useful for
development.

First, enter the virtual environment:

.. code:: shell

    workon voseq_environment

Then, create tables for the database:

.. code:: shell

    cd /path/to/Voseq/
    make migrations

And import test data for your database:

.. code:: shell

    make test_import

Start a test server
===================

In Linux start elasticsearch as a service, then enter the virtual environment and then start the server:

.. code:: shell

    sudo service elasticsearch start
    workon voseq_environment
    cd /path/to/Voseq
    make serve

In macOSX if you do not have the ``service`` command, run
``elasticsearch`` in the background and then start the server (\*):

.. code:: shell

    elasticsearch -d
    cd /path/to/Voseq
    make serve

\* *Note that if you did not check to Start Postgres automatically after
login, you first have to go to Applications and start it manually from
there by clicking on the Postgres.app. Do this before running the
server.*

You now have a local webserver running. You can access it by opening this URL in your web browser: ``http://127.0.0.1:8000/`` and try all the buttons to see if they all work! Also notice the debug bar on the right of the screen where you can check if all the configurations are correct.

Migrate VoSeq database
======================

If you have an existing Voseq 1.x database and want to migrate, you need to dump your MySQL database into a XML file:

.. code:: shell

    cd /path/to/Voseq/
    mysqldump --xml voseq_database > dump.xml

Then use our script to migrate all your VoSeq data into a PostGreSQL
database.

.. code:: shell

    make migrations
    python voseq/manage.py migrate_db --dumpfile=dump.xml --settings=voseq.settings.local

If you have used a prefix for your tables in the old VoSeq, you can optionally input this as an
argument for the import script:

.. code:: shell

    python voseq/manage.py migrate_db --dumpfile=dump.xml --prefix=voseq_ --settings=voseq.settings.local


It might issue a warning message:

::

    WARNING:: Could not parse dateCreation properly.
    WARNING:: Using empty as date for `time_edited` for code Your_Vocher_Code

It means that the creation time for your voucher was probably empty or
similar to ``0000-00-00``. In that case the date of creation for your
voucher will be empty. This will not cause any trouble when running
VoSeq. You can safely ignore this message.

Create an index for all the data in your database:

.. code:: shell

    make index

Set-up a publicly available web server
======================================

To make VoSeq available to multiple users, you will have to set-up a publicly available web server. There are several options to do this, for example using nginx and gunicorn (best performance) or Apache and WSGI (more suitable for hosting multiple websites).

Instructions for how to do this will follow later, but the DigitalOcean tutorials may be of use for now:

`Apache and WSGI <https://www.digitalocean.com/community/tutorials/how-to-run-django-with-mod_wsgi-and-apache-with-a-virtualenv-python-environment-on-a-debian-vps>`__

`Nginx and Gunicorn <https://www.digitalocean.com/community/tutorials/how-to-install-and-configure-django-with-postgres-nginx-and-gunicorn>`__

Administrate the server
=======================

Optionally if you want to add items/vouchers to your database
interactively, you need to create an administration account. Run the
following command and provide the requested information:

.. code:: shell

    make admin


Some features of VoSeq need to be run periodically. You can setup cronjobs to execute some commands once a day or every 2 hours depending on your needs:

* Update the database index for the simple and advanced search functions:

.. code:: shell

    python voseq/manage.py update_index --settings=voseq.settings.local

* Update some voucher and gene statistics for your installation of VoSeq:

.. code:: shell

    make stats


.. |Waffle| image:: https://badge.waffle.io/carlosp420/voseq.png?label=in%20progress&title=In%20Progress
   :target: https://waffle.io/carlosp420/voseq
   :alt: 'Stories in Progress'
.. |Dependency Status| image:: https://gemnasium.com/carlosp420/VoSeq.svg
   :target: https://gemnasium.com/carlosp420/VoSeq
.. |Coverage Status| image:: https://img.shields.io/coveralls/carlosp420/VoSeq.svg
   :target: https://coveralls.io/r/carlosp420/VoSeq?branch=master
.. |Build status| image:: https://ci.appveyor.com/api/projects/status/0ba440vjw8811845/branch/master?svg=true
   :target: https://ci.appveyor.com/project/carlosp420/voseq/branch/master
.. |Build Status| image:: https://travis-ci.org/carlosp420/VoSeq.svg
   :target: https://travis-ci.org/carlosp420/VoSeq
.. |Landscape| image:: https://landscape.io/github/carlosp420/VoSeq/master/landscape.svg
   :target: https://landscape.io/github/carlosp420/VoSeq/master
   :alt: Code Health
.. |Docs| image:: https://readthedocs.org/projects/voseq/badge/?version=latest
   :target: http://voseq.readthedocs.org/en/latest/
   :alt: Documentation Status


.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/carlosp420/VoSeq
   :target: https://gitter.im/carlosp420/VoSeq?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge