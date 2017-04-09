^^^^^^^^^^^^^^^^^^^
How to update Voseq
^^^^^^^^^^^^^^^^^^^

Versions 2.0.0 and later
""""""""""""""""""""""""
The easiest way to update VoSeq (that does not require new install of software or database):

* `Download the new files from github <https://github.com/carlosp420/VoSeq/tags>`_.

    * Unpack the new files to your **webserver directory** (htdocs, webserver, etc).
    * Rename your old "in use" VoSeq folder something like, ``VoSeq_old`` or similar (e.g. ``VoSeq`` -> ``VoSeq_old``).
    * Give the newly downloaded VoSeq folder the same name as the old one had (e.g. ``VoSeq``).
    * Copy the file ``conf.json`` (in main folder) from the old version (e.g. ``VoSeq_old``) to the new version (e.g. ``VoSeq``).


Versions 1.7.4 and earlier
""""""""""""""""""""""""""

As VoSeq has been rewritten completely you will need to download the source code
and install it from scratch (See :ref:`intro-install`).

If you have been using versions 1.7.4 and earlier you need to migrate your data
from MySQL to postgreSQL (See :ref:`migrate-mysql`).


