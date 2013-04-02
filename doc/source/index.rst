.. VoSeq documentation master file, created by
   sphinx-quickstart on Mon Apr  1 22:15:34 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to VoSeq's documentation!
=================================

Contents:

.. toctree::
   :maxdepth: 2


..
	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><
	Page 1.

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





..
	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><
	Page 2.

---------------
Getting Started
---------------
Once you have successfully downloaded VoSeq, you can find out how to:

* :ref:`install_in_linux`
* [[Install in Mac OS X]]
* [[Install in Windows 7 / Vista / XP]]
* [[Quick Guide]] to get started with VoSeq.


.. _install_in_linux:

^^^^^^^^^^^^^^^^
Install in Linux
^^^^^^^^^^^^^^^^
Before installing VoSeq, you need to install in your computer a web server 
(such as `Apache <http://httpd.apache.org/>`_) and the relational database
`MySQL <http://www.mysql.com/>`_.


"""""""""""""""""
Required software
"""""""""""""""""

* Web server with PHP 5.0 or higher (http://www.php.net/manual/install.unix.php). **Compile it with the library CURL**, which is needed to do BLASTs against GenBank.

	* Apache HTTP Server
	* PHP
* A MySQL server 5.0 or higher (see http://www.mysql.com)
* GD library

These instructions assume that your are using Linux and Apache, and have installed
`LAMP <http://en.wikipedia.org/wiki/LAMP>`_ (Linux, Apache, MySQL and PHP on your
computer).

#. Compile PHP with support for the graphics library GD. More info `here <http://www.php.net/manual/en/image.installation.php>`_.
#. Download VoSeq: `Download from github <https://github.com/carlosp420/VoSeq/tags>`_.
#. Unzip the source files in some directory: ``unzip VoSeq_X.Y.Z.zip``
#. If you are not a Linux Guru and you have `WinRAR <http://www.rarlab.com/>`_ (like WinZip but works with gzipped files) on your Windows system you can cheat a little bit here. You can download the file to your Windows machine, use WinRAR to unzip the gzipped file into a directory in Windows and then use an FTP program like `WinSCP <http://winscp.net/eng/index.php>`_ to transfer the entire VoSeq directory for you to a commercial server for example.
#. Move the directory into your web directory: e.g. ``mv VoSeq /usr/local/apache2/htdocs/myVoSeq`` or ``mv VoSeq public_html/myVoSeq`` or use your FTP software to do this for you.
#. To run the installation script, you'll need to temporarily make your myVoSeq directory writable by the web server. The simplest way to do this on a Unix/Linux system is to make it world-writable by typing: ``chmod -R 777 myVoSeq``. To do this into a commercial server you will need a telnet client like `PuTTY <http://www.chiark.greenend.org.uk/~sgtatham/putty/>`_ on your system.
#. At this point you should have Apache and MySQL running (this varies between distributions and setups, see their documentations for details).
#. Go to your web browser and surf into the VoSeq installation directory (under ``htdocs`` or ``public_html`` folders of Apache). It will direct you to the config script (if it doesn't, just load up the ``http://localhost/myVoSeq/index.php`` file. Fill out the forms.
#. If all goes well, the installer will create a configuration file named ``conf.php`` in your myVoSeq installation directory. This file will contain all the important variables and information needed to run VoSeq in your system. 


.. _blast-plugin:

^^^^^^^^^^^^
BLAST plugin
^^^^^^^^^^^^
against all others in your VoSeq database (ssection for
xxxxxxxxx


