.. _development:

^^^^^^^^^^^^^^^^^^^^
Development of VoSeq
^^^^^^^^^^^^^^^^^^^^

VoSeq has been constructed using the framework for web applications Django. In
this document, you will find some guidelines about the structure of the software
of VoSeq.

----------------
Dataset Creation
----------------

Module: ``voseq.create_dataset``

This module uses the Python library Dataset-creator_ which uses some specific
parameters to create datasets. These parameter values should be used in VoSeq,
starting from the forms and models behind the graphical interface as well as
the internal modules that handle the dataset creation.

- partitioning (str):
    Partitioning scheme:  ``by gene`` (default), ``by codon position`` (each) and ``1st-2nd, 3rd``.


.. _Dataset-creator: https://github.com/carlosp420/dataset-creator
