.. VoSeq documentation master file, created by
   sphinx-quickstart on Mon Apr  1 22:15:34 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#########################################
Welcome to VoSeq |version| documentation!
#########################################

============
Getting help
============

* Try the :doc:`FAQ <faq>` -- with answers to common questions.
* If you need help regarding installation or usage of th application, please
  contact `Carlos Peña <mycalesis@gmail.com>`_.
* You can also subscribe to VoSeq's discussion list on
  `Google Groups <https://groups.google.com/d/forum/voseq-discussion-list>`_.

.. toctree::
   :maxdepth: 2
   :hidden:

   news
   faq


Welcome to **VoSeq**, a database to store voucher and DNA sequence data for
phylogenetic analysis. It's a unique database that enables to digitize
biological data of museum specimens and molecular data such as DNA sequences,
primers and genes.

VoSeq has tools that facilitate the batch upload of lots of voucher data and
DNA sequences with a few clicks. It also has
`BLAST <http://en.wikipedia.org/wiki/BLAST/>`_ capabilities, meaning that you
can find out whether one particular
DNA sequence is most similar to other sequence in `NCBI GenBank
<http://www.ncbi.nlm.nih.gov/genbank/>`_. You can also BLAST your sequence
against all others in your VoSeq database (see :ref:`blast-plugin` section for
details).

VoSeq is written in `Python <https://www.python.org/>`_. It uses
`PostgreSQL <http://www.postgresql.org/>`_ as database back-end and it is
designed to run either locally in your own computer or on a remote (commercial)
server service.

==========
Start here
==========

.. toctree::
   :hidden:

   intro/overview
   intro/install

   usage/adding-vouchers
   usage/adding-genes
   usage/adding-sequences
   usage/create-taxonset
   usage/create-datasets
   usage/search
   usage/upload-voucher-photos
   usage/create-excel-table
   usage/update-voucher

   update-voseq
   plugins

:doc:`intro/overview`
   Find out what VoSeq can do. It might be right for you.

:doc:`intro/install`
   Install VoSeq on your computer or server.

:ref:`configure`
   Create a ``conf.json`` file to specify your settings.


.. _citing-voseq:
=================
How to cite VoSeq
=================
If you think VoSeq is useful and you happen to use it during your work, it
would be great if you cite us as a source:

.. note::

    Peña, C. & Malm, T. **2012**. VoSeq: a Voucher and DNA Sequence Web Application.
    *PLOS ONE*, 7(6): e39071.
    `doi: 10.1371/journal.pone.0039071 <http://dx.doi.org/10.1371/journal.pone.0039071>`_

===========
Using VoSeq
===========

:doc:`usage/adding-vouchers`
    Learn how to create records for voucher specimens in VoSeq.

:doc:`usage/adding-genes`
    Learn how to add data about sequenced genes in VoSeq.

:doc:`usage/adding-sequences`
    Learn how to upload DNA sequences to VoSeq.

:doc:`usage/create-taxonset`
    Learn how to group set of vouchers in TaxonSets.

:doc:`usage/create-datasets`
    Learn how to create dataset files for phylogenetic software.

:doc:`usage/search`
    Learn how to do simple and advanced searches.

:doc:`usage/upload-voucher-photos`
    Learn how to upload voucher pictures in VoSeq.

:doc:`usage/create-excel-table`
    Learn how to export data ready for Excel tables.

:doc:`usage/update-voucher`
    Learn how to update voucher and sequence data.

===================
Management of VoSeq
===================


:doc:`update-voseq`
    Upgrading VoSeq from versions 1.7.4 and greater can be done rather simply if you
    downloaded VoSeq using ``git``.

:doc:`plugins`
    Some tools in VoSeq can be enabled by setting up plugins.
