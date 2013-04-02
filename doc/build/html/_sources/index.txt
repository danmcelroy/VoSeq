.. VoSeq documentation master file, created by
   sphinx-quickstart on Mon Apr  1 22:15:34 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to VoSeq's documentation!
=================================

Contents:

.. toctree::
   :maxdepth: 2


-----------
Hello there
-----------

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

VoSeq is written mainly in `PHP <http://www.php.net/>`_. It uses 
`MySQL <http://www.mysql.com>`_ as back-end and it is designed to run in a
local server (for example by installing `Apache <http://httpd.apache.org/>`_
on your computer) or to run on any commercial server service.

.. image:: images/intro1.png
   :align: center
   :width: 240px
.. image:: images/create_taxonset2_small.png
   :align: center
.. image:: images/create_dataset_small.png
   :align: center


^^^^^^^^^^^^^^^^^
How to cite VoSeq
^^^^^^^^^^^^^^^^^
If you think VoSeq is useful and you happen to use it during your work, it
would be great if you cite us as a source:

* Peña, C. & Malm, T. **2012**. VoSeq: a Voucher and DNA Sequence Web Application. *PLOS ONE*, 7(6): e39071.  `doi <http://dx.doi.org/10.1371/journal.pone.0039071>`_

^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Help and contact information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you need help regarding installation or usage of th application, please
contact `Carlos Peña <mycalesis@gmail.com>`_ or `Tobias Malm <tobemalm@gmail.com>`_.

You can also subscribe to VoSeq's discussion list on `Google Groups <https://groups.google.com/d/forum/voseq-discussion-list>`_.


.. _blast-plugin:
------------
BLAST plugin
------------
xxxxxxxxx

