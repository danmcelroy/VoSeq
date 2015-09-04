.. _migrate-mysql:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Migrate VoSeq's MySQL database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VoSeq's versions 1.7.4 and earlier used the MySQL database to store data.
Starting with version 2.0.0 the database to be used is postgreSQL.

Migrating the data can be easily done from the command line:

* Step 1. Export your data in MySQL as an XML formatted file:

.. code-block:: shell

    mysqldump --xml voseq_database -u user > dump.xml

* Step 2. After you have downloaded the source code of VoSeq (versions 2.0.0 and later)
  and configured your postgreSQL installation (See :ref:`intro-install`), you can
  import your data with the following commands:

.. code-block:: shell

    cd /path/to/VoSeq/
    make migrations
    python voseq/manage.py migrate_db --dumpfile=/path/to/dump.xml --settings=voseq.settings.local

If you have used a prefix for your tables in the old VoSeq, you can optionally
input this as an argument for the import script:

.. code-block:: shell

    python voseq/manage.py migrate_db --dumpfile=dump.xml --prefix=voseq_ --settings=voseq.settings.local
