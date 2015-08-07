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

:doc:`intro/overview`
   Find out what VoSeq can do. It might be right for you.

.. image:: images/intro1.png
   :align: center
   :width: 240px

.. image:: images/create_taxonset2_small.png
   :align: center
   :width: 357px

.. image:: images/create_dataset_small.png
   :align: center
   :width: 329px

.. _MySQL: http://www.mysql.com

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
* :ref:`install_in_mac`
* :ref:`install_in_windows`
* :ref:`quick_guide` to get started with VoSeq.


.. _install_in_linux:

^^^^^^^^^^^^^^^^
Install in Linux
^^^^^^^^^^^^^^^^
Before installing VoSeq, you need to install in your computer a web server
(such as `Apache <http://httpd.apache.org/>`_) and the relational database
`MySQL`_.


"""""""""""""""""
Required software
"""""""""""""""""

* Web server with PHP 5.0 or higher (http://www.php.net/manual/install.unix.php). **Compile it with the library CURL**, which is needed to do BLASTs against GenBank.

	* Apache HTTP Server
	* PHP
* A MySQL server 5.0 or higher (see http://www.mysql.com)
* GD library

.. note:: These instructions assume that your are using Linux and Apache, and have installed `LAMP <http://en.wikipedia.org/wiki/LAMP>`_ (Linux, Apache, MySQL and PHP on your computer).

#. Compile PHP with support for the graphics library GD. More info `here <http://www.php.net/manual/en/image.installation.php>`_.
#. Download VoSeq: `Download from github <https://github.com/carlosp420/VoSeq/tags>`_.
#. Unzip the source files in some directory: ``unzip VoSeq_X.Y.Z.zip``
#. If you are not a Linux Guru and you have `WinRAR <http://www.rarlab.com/>`_ (like WinZip but works with gzipped files) on your Windows system you can cheat a little bit here. You can download the file to your Windows machine, use WinRAR to unzip the gzipped file into a directory in Windows and then use an FTP program like `WinSCP <http://winscp.net/eng/index.php>`_ to transfer the entire VoSeq directory for you to a commercial server for example.
#. Move the directory into your web directory: e.g. ``mv VoSeq /usr/local/apache2/htdocs/myVoSeq`` or ``mv VoSeq public_html/myVoSeq`` or use your FTP software to do this for you.
#. To run the installation script, you'll need to temporarily make your myVoSeq directory writable by the web server. The simplest way to do this on a Unix/Linux system is to make it world-writable by typing: ``chmod 777 myVoSeq``. To do this into a commercial server you will need a telnet client like `PuTTY <http://www.chiark.greenend.org.uk/~sgtatham/putty/>`_ on your system.
#. At this point you should have Apache and MySQL running (this varies between distributions and setups, see their documentations for details).
#. Go to your web browser and surf into the VoSeq installation directory (under ``htdocs`` or ``public_html`` folders of Apache). It will direct you to the config script (if it doesn't, just load up the ``http://localhost/myVoSeq/index.php`` file. Fill out the forms.
#. If all goes well, the installer will create a configuration file named ``conf.php`` in your myVoSeq installation directory. This file will contain all the important variables and information needed to run VoSeq in your system.


.. _install_in_mac:

^^^^^^^^^^^^^^^^^^^
Install in Mac OS X
^^^^^^^^^^^^^^^^^^^
We have successfully installed VoSeq in a MacBook OS X Lion. It appears that the Mac
operative systems **come already with Apache and PHP installed**. However you will
need to enable Apache to read and run PHP files.


"""""""""""""""""""""""""""""""""""""""""""""""""""""
To connect Apache and PHP so that they work together:
"""""""""""""""""""""""""""""""""""""""""""""""""""""

#. Edit Apache's configuration text file:

    * ``sudo nano /etc/apache2/httpd.conf``
#. Make sure that the line: ``LoadModule php5_module     libexec/apache2/libphp5.so``  is in the file and it is not commented (there is no # symbol at the beginning of the line).
#. Find the section ``<IfModule mime_module>`` and write the following line ``AddType application/x-httpd-php .php`` so that Apache will run any file with the extension .php as a script and will not show it as plain text.


"""""""""""""
Install MySQL
"""""""""""""
Unfortunately Mac OS X systems don't come with MySQL installed. You can download it from here:

#. Download MySQL from here: http://dev.mysql.com/downloads/mysql/5.1.html

    * Download the ``.dmg`` package according to your systems specifications (32 bits or 64 bits).
#. You might also want to install MySQL GUI Tools http://dev.mysql.com/downloads/gui-tools/5.0.html
#. The following is a quick guide to installling MySQL on your computer. **It is not comprehensive and you will find much more info in the documentation for installing MySQL here**: http://dev.mysql.com/doc/mysql-macosx-excerpt/5.5/en/index.html
#. Unpack and install both pieces of software. Make sure you install the package, in my case, ``mysql-5.1.60-osx10.6-x86_64.pkg`` and ``MySQLStartupitem.pkg``
#. Start the MySQL server by typing in the terminal: ``sudo /Library/StartupItems/MySQLCOM/MySQLCOM start``
#. Create a password for the user **root** by typing: ``/usr/local/mysql/bin/mysqladmin -uroot password 'myownpassword'``


"""""""""""""
Install VoSeq
"""""""""""""
#. To start Apache, go to System Preferences>Sharing> and tick Web Sharing to start your web server. Your assigned folder to host your webpages and VoSeq installation is the folder Sites in your Home directory: ``/Users/YourName/Sites``. You will need to place there the source files of ``VoSeq_X.Y.Z.zip``
#. You need to click the button "create personal share folder" to create the folder "Sites".
#. Open a Terminal: go to Applications>Utilities>Terminal. In the Terminal window, type ``cd ~/Sites`` to go to the folder where the file ``Voseq_X.Y.Z.zip`` should be.
#. Unpack the contents by typing ``unzip VoSeq_X.Y.Z.zip``
#. Start the MySQL server: ``sudo /Library/StartupItems/MySQLCOM/MySQLCOM start``
#. Go to your web browser and point it to the VoSeq installation directory: ``http://localhost/~YourName/VoSeq``. It will direct you to the config script. Fill out the forms.
#. If all goes well, the installer will create a configuration file named ``conf.php`` in your VoSeq installation directory. This file will contain all the important variables and information needed to run VoSeq in your system.
#. If during installation, VoSeq cannot connect to MySQL server, you might need to modify your ``/usr/local/mysql/support-files/my-large.cnf`` file parameters:

    * Modify the lines ``/var/mysqld/mysqld.sock`` to  this ``/tmp/mysql.sock``
    * Save the file as ``/etc/my.cnf``



.. _install_in_windows:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Manual Install in Windows 7 / Vista / XP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow these instructions to install Apache, PHP and MySQL and lastly VoSeq on Windows 7 / Vista / XP systems - its not as hard as it looks!

""""""
Apache
""""""

#. **Download and install "Apache2.x"** (tested on 2.2.x) http://httpd.apache.org/ as recommended, preferrably use "localhost" as Network Domain and Server Name. Start the service and try it out by opening http://localhost in your web browser - the output should be **It works!**.

.. note:: Notice that Apache will want to use port ``0.0.0.0:80``, which may be used by other programs, if Apache doesnt start (may say something about port occupied), write ``netstat -nab`` in Terminal and check if some other process is using that adress - then close that process if appropriate.

"""
PHP
"""

#. **Download "PHP 5.x.zip"** (tested on version 5.2.17) http://windows.php.net/download/. We recommend that you download the  ``VC6 Thread Safe`` version if using Apache. Unpack the ``PHP5.x.zip`` file to a folder named PHP (ex. ``C:\PHP`` or ``C:\Program Files\PHP``). Then copy the ``php.ini-???`` to ``C:\WINDOWS`` and rename it ``php.ini``. (``???`` can be dist, production or development).
#. Open the apache configuration file ``httpd.conf`` in a text editor (found in the ``C:\Program Files\Apache Software Foundation\Apache2.2\conf`` folder after standard install).
#. Add the following 4 lines at the end of the ``LoadModule`` section (now assuming php installed to ``C:\PHP`` otherwise change this to correct installation folder)::

    LoadModule php5_module "c:/PHP/php5apache2_2.dll"
    AddHandler application/x-httpd-php .php
    # configure the path to php.ini
    PHPIniDir "c:/windows"

#. Add a file called ``info.php`` containing ``<?php phpinfo();?>`` to the ``C:\Program Files\Apache Software Foundation\Apache2.2\htdocs`` folder.
#. Restart your Apache Server to confirm changes: "Start > All Programs > Apache HTTP Server 4.2.4 > Control Apache Server > Restart".
#. Open up your web browser and type in: http://localhost/info.php. If you get a page with blue tables containing PHP and Apache info, then **installation is successful!**
#. Finish installing PHP by modifying your PHP Configuration File (``C:\WINDOWS\php.ini``) in a text editor:

    * Find the line containing: (Delete the "``;``" at the beginning of the lines)

        * ``;extension_dir = "./"`` and change it to
        * ``extension_dir = "C:\php\ext"``

    * and the line containing:

        * ``;session.save_path = "/tmp"``" and change it to
        * ``session.save_path = "C:\WINDOWS\temp"``


""""""""""""""""""""""""
Enable the curl protocol
""""""""""""""""""""""""

Curl is needed to get the Flickr plugin to work and enable VoSeq to interact with other databases.

#. Copy the file ``php_curl.dll`` from the folder ``C:\PHP\ext`` into the folder ``C:\WINDOWS\system32``
#. Remove the semicolon ``;`` from the line ``;extension=php_curl.dll`` in your file ``C:\WINDOWS\php.ini``
#. Restart the apache server.


"""""
MySQL
"""""

#. **Download and install MySQL** (tested on 5.5) from http://dev.mysql.com/downloads/mysql/ with typical install - check the "skip Sign-Up" and '"Configure the MySQL server now" boxes when they arrive. Finish installation.
#. The MySQL Server Instance Configuration Wizard should appear.

    * Click "next" ->
    * Select "Detailed Configuraton" and click "next" ->
    * Select "Developer Machine" and click "next" ->
    * Select "Multifunctional Database" and click "next" -> click "next" ->
    * Select "Decision support (DSS)/OLAP" and click "next" ->
    * Check "Enable TCP/IP Networking"
    * "Port Number" should be set to "3306" and
    * Check "Enable strict mode", click "next" ->
    * Select "Standard Character Set" and click "next" ->
    * Check "Install As Windows Service, set the name to "MySQL" and check "Launch the MySQL Server automatically
    * Make sure that the "Include Bin Directory in Windows Path" **is NOT checked**.
    * Click "Next". -> Check the box that says "Modify Security Settings".
    * Enter a password for the default "root" account, and confirm the password in the box below.
    * **Do NOT check the boxes** "Enable root access from remote machines" or "Create An Anonymous Account".
    * Click "Next" -> Click "Execute" and let it finish.
    * Click "Finish". Now MySQL should be installed.

#. **To enable PHP to use the MySQL databases**, open the ``php.ini`` (``C:/WINDOWS/php.ini``) file in your text editor and find the line ``;extension=php_mysql.dll``. Delete the "``;``" at the beginning of the line and save the file.
#. Add the PHP directory to Windows PATH - To do this, click:

    * Start > My Computer > Properties > Advanced > Environment Variables.
    * Under the second list (System Variables), there will be a variable called "Path".
    * Select it and click "Edit". Add "``;C:\php``" (or your own path to PHP if installed as other) to the very end of the string and click "OK".

#. Restart your computer and try out the database.
#. (Optional) In order to easily makes changes or additions in your database download and try out the `MySQL Workbench <http://dev.mysql.com/downloads/workbench/5.2.html>`_



"""""""""""""
Install VoSeq
"""""""""""""

#. Download and unzip the file ``Voseq_VersionNumber.zip`` in the Apache folder (rename the new folder if necessary):

    * ``C:\Program Files\Apache Software Foundation\Apache2.2\htdocs``

#. Point your web browser to the address (that is - localhost + the name of your VoSeq folder): ``http://localhost/VoSeq_VersionNumber`` and follow the instructions for installing the software.
#. If all goes well, the installer will create a configuration file named ``conf.php`` in your VoSeq installation directory. This file will contain all the important variables and information needed to run VoSeq in your system.



.. _install_with_xampp:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Install in Windows with XAMPP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you dont want to install Apache, MySQL and PHP manually you may want to try using a 3-rd party installer.
XAMPP installs all three as well as the extra protocols needed for PHP.

"""""""""""""
Install XAMPP
"""""""""""""

You can download XAMPP installer at http://www.apachefriends.org/en/xampp.html
Then install it with the installer (tested with version 1.8.1).
With XAMPP MySQL is installed without password, **for security you need to create a new password for MySQL** (as well as for the XAMPP web directory which by default is accessible for everyone that know your IP adress, though you may still be somewhat protected behind a router (`read here <http://www.apachefriends.org/en/xampp-windows.html#1221>`_).

* Goto ``localhost/security`` and check your security level and set passwords!

"""""""""""""
Install VoSeq
"""""""""""""

#. Download and unzip the file ``Voseq_VersionNumber.zip`` in the XAMPP/htdocs directory (rename the new folder if necessary):

    * ``C:\XAMPP\htdocs``

#. Point your web browser to the address (that is - localhost + the name of your VoSeq folder): ``http://localhost/VoSeq_VersionNumber`` and follow the instructions for installing the software.
#. If all goes well, the installer will create a configuration file named ``conf.php`` in your myVoSeq installation directory. This file will contain all the important variables and information needed to run VoSeq in your system.

Configuration files after XAMPP install can be seen :ref:`xampp_config`.

.. note::
    * This has only been tested quickly, and may not work for all computer systems!

        * We welcome all feedback for this type of installation!

    * If you already have MySQL install XAMPP SHOULD not overwrite your existing databases, but precaution is a virtue (or something...) and we advice making backups of stored data before installation. (We can not be held responsible for any loss of data)
    * More information regarding XAMPP for windows is found here: http://www.apachefriends.org/en/xampp-windows.html


.. _xampp_config:

""""""""""""
XAMPP config
""""""""""""

.. image:: images/XAMPP_configs.png
   :align: center
   :width: 571px

..
	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><
	Page 3.

.. _quick_guide:





.. _adding_vouchers:

^^^^^^^^^^^^^^^
Adding vouchers
^^^^^^^^^^^^^^^

After successful installation, the first thing to do is to add records (vouchers, or specimens). You can add a single record by going to the **Administrator interface** and clicking on the link **Add new record**.

.. image:: images/add_new_record.png
   :align: center
   :width: 446px

The most important information to enter is the **code** of voucher, which has to be unique. VoSeq will refuse to accept duplicate codes and will issue error message if this happens. Another necessary field is the **Genus** entry, while all other fields are optional.

You can also upload a batch of records using the tool **Upload batch sequences/vouchers**. You will be shown a page to batch-upload sequences. By clicking the button **Upload vouchers instead** you will see instructions on how to upload specimen data. You can quickly import voucher data from a table in MS Excel by copying and pasting into the text area, provided that you use the right field headers.

.. image:: images/batch_record_upload.png
   :align: center
   :width: 902px



.. _adding_genes:

^^^^^^^^^^^^
Adding genes
^^^^^^^^^^^^

The second thing after adding vouchers that should be done is to create **new genes** or "alignments" for your database. This must be done in order to be able to **add sequences** to the database.
You can add a single record by going to the **Administrator interface** and clicking on the link **Add/edit/view gene information**, followed by **Add gene**'.

.. image:: images/add-edit_gene.png
   :align: center
   :width: 540px

The most important and the only field that is obligatory is the **gene code** field, this will be the name of your gene when using the database. This could be a simple short version (e.g. **COI** for **Cytochrome oxidase I**) or any other name (no spaces allowed, but **_** are ok).
You can for example create a gene code for aligned data, say the barcode version of COI of 658 bp could be named COI_658 or similar. Other genecodes could be made for unaligned sequences (e.g. COI_raw).

First the gene/alignment have to be **specified if aligned or not**. For example may raw sequences be set to **no** (and these may be retrieved as FASTA-files, whereas if you want to build other datasets (Nexus, PHYLIP, TNT) they need to be set to aligned. If set to **no**, then other information regarding reading frames and such will be ignored!

Then (if aligned) you should include the **length** of an aligned gene - this will be used for dataset creation and will there warn for sequences longer than the specified length.

You may also add a **description** for the gene - these should be the **full name of the gene** (e.g. **Cytochrome oxidase I**) - as this field is used for example in creating tables to submit to **GenBank**.

Aligned genes may be set as **protein-coding** for **additional prot-coding gene functionality** when retrieveing datasets for example (e.g. position choices, translation).

For aligned protein coding (=yes) genes you may choose to specify the **reading frame** as well as the genetic code for translation - this will be used for dataset creation and is a must f you want to partition your genes according to **codon positions** or **amino acids**.

**Introns** may be added - enter number of introns in your alignment and click 'update introns', that will give new fields for entering starting and finsihing positions for your introns. (Remember that positions in an alignment here is counted from 1 and upwards).


.. _adding_sequences:

^^^^^^^^^^^^^^^^
Adding sequences
^^^^^^^^^^^^^^^^

In the **Administrator interface**, the tool **Upload batch sequences/vouchers** allows you to upload DNA sequences into VoSeq. Along with the DNA sequences, you have to upload the required fields **gene code** and **voucher code**, optionalyl the primer names, laborator and creation date. Each sequence and its related data goes into one line, with fields separated by tabs. If you have your data in a spreadsheet such as MS Excel, you can copy and paste the data into the text area.

* It is important to use the same headers provided in the text area.
* It is also necessary that that the **code** of each sequence matches the **code** of voucher specimens that had been uploaded into VoSeq. This is the unique identifier that is used to connect the voucher data and their sequences.

(Aligned sequences should for best use of the database have missing data coded as questionmark (?) and gaps as a dash (-))

.. image:: images/batch_seq_upload.png
   :align: center
   :width: 902px



.. _`create_taxonset`:

^^^^^^^^^^^^^^^
Create taxonset
^^^^^^^^^^^^^^^

Taxonsets is a way to make a list of taxa that are being used for a specific project or analysis. A Taxonset is just a list of voucher codes. By having Taxonsets, you can quickly create datasets and tables for them.

If you have not set Taxonsets you will need to type specimen codes everytime you create a dataset. Instead, if you have a Taxonset for a particular project, you could select it when creating Tables for manuscripts.

Go to the **Administrator interface** and click on the link **Add/edit/view Taxon sets**.

A taxonset must have a name in order to be saved and usable later!

You can create a **Taxonset** by entering a list of specimen codes, each separated by a return:

.. image:: images/create_taxonset1.png
   :align: center
   :width: 518px

Or by browsing the data in VoSeq and choosing the specimens you are interested in by marking them in the **X** field:

.. image:: images/create_taxonset2.png
   :align: center
   :width: 792px

Here you can sort the table according to choosen information (taxonomic level, code, X-marked or not), as well as choose genes to display information of.
If you have choosen one or several genes, you can sort the table according to sequence availability for selected genes.
You can also press **mark all** or **unmark all** to add or remove X's to or from each taxa that are displayed (works well with filtering).
In order to perform a filtering or after selection of a new genecode you must press **Sort/Filter** to proceed. Your already marked taxa will be remembered.

After completing your selection of taxa and adding name and descriptions - press the **Add dataset** button to save it.
If you are updating an already existing taxonset - press **Update taxon set**.



.. _create_datasets:

^^^^^^^^^^^^^^^
Create datasets
^^^^^^^^^^^^^^^

We believe that one of VoSeq's important features is the **capability to create dataset files of molecular sequences that are ready-to-run in phylogenetic software** such as **MrBayes, TNT, PAUP, RaXML**, etc.

Now that you have voucher and sequence data in your installation of VoSeq, you might want to create datasets for analysis of sequences in phylogenetic software.

In the **user interface**, you will find under the **Tools** section the link **Create new dataset**. You will be shown a page to select the sequences you want by entering the **voucher codes** and **gene codes**. You can select your data to be in several formats (FASTA, NEXUS, etc), choose between codon positions, as well as choosing what information your taxon names should include.

This will create a **ready-to-use** data set for analyses!

.. image:: images/create_dataset_new.png
   :align: center
   :width: 534px

The **Outgroup** field, if needed, should include the voucher code for the chosen outgroup taxa.

Codon positions
	Marking **1st**, **2nd** or **3rd** and unmarking **all** positions will create a dataset with only the chosen position(s) for all genes.
	**Special** will take you to a new page where you will be able to choose which codon positions to include for each gene.
	Marking **amino acids** will tranlate **protein-coding** genes with a set **genetic code** , the others will be treated as normal dna, i.e. making "mixed" datatype in Nexus for MrBayes, and setting partitions correctly in Nexus and PHYLIP format.
	Note that codon position choices as well as translation to aminoacids are only able to function if the user have specified a **reading frame** for the chosen gene(s) (see :ref:`adding_genes`).

Partition by (position)
	Here you can choose how to do partitioning for each gene.
	**as one** will create one partition per gene, regardless of which codon positions you include.
	**each** will create a partition per codon position, whereas **1st-2nd, 3rd** will create one combined partition for the 1st and 2nd positions and one separate for the 3rd codon positions.
	Note that **each** and **all** are only possible to process with a per gene specified **reading frame** (see :ref:`adding_genes`).

You can also chose to **omit taxa from a taxonset that contains less than a specific number sequences**. Say you have a 10 gene data set and want to remove all taxa with 5 or less! Easy! Just eneter a minimum number of genes!

If you have introns in your alignment you can choose to include or remove them from the output data set. If included they will be treated as separate data blocks and partitions for the Nexus and PHYLIP outputs!

The voucher codes can be entered one by one (separated by return) in the text area or you could create a :ref:`create_taxonset` (a list of voucher codes for a specific project).

.. note:: As of version 1.5.0, protein-coding ability, aligned or not, introns and genetic code will be set for each gene/alignment in the admin gene section!



.. _my_search:

^^^^^^
Search
^^^^^^

You can search for records by queries using single fields or any combination of them. The autocomplete dropboxes will help you query existing data easily.
This can be done in both the **user interface** and the **administrator interface** - where the latter have more options to search (e.g. record history).

.. image:: images/search.png
   :align: center
   :width: 590px



.. _upload_voucher_photos:

^^^^^^^^^^^^^^^^^^^^^
Upload voucher photos
^^^^^^^^^^^^^^^^^^^^^

In the **Administrator interface** you will see that some records have the link **Picture missing**. By clicking on this link, you will be able to upload a photo for that voucher.

.. image:: images/picture_missing.png
   :align: center
   :width: 376px

If you want to replace an existing picture with another, you will need to click the **Change picture** icon.

After you upload your photo, VoSeq will automatically **post the picture in Flickr** and save the necessary URL addresses in the MySQL database. Thus, you will see your photo in the corresponding voucher page.

If you have not enabled the :ref:`flickr_plugin`, VoSeq will instruct you how to do this.

**If you don't want to use Flickr**, you can host your photos locally on your own server or computer. For this you will need to edit a line in your ``conf.php`` file:

* Change the line:

	* ``$photos_repository = 'flickr';`` to this one:
	* ``$photos_repository = 'local';``

Starting with version 1.5.0, VoSeq can host many photos for each voucher. Photos can be added in the voucher page using the administrator interface of VoSeq. You can delete photos individually by clicking on it's "trash" icon.

.. note:: If you have more than two photos for voucher, all additional photos will appear at the bottom of the voucher page (see image below).

.. image:: images/additional_voucher_photos.png
   :align: center
   :width: 433px

.. _create_excel_table:

^^^^^^^^^^^^^^^^^^
Create Excel table
^^^^^^^^^^^^^^^^^^
You can create a MS Excel table with specimen codes, genus and species names, genes used in analysis along with their accession numbers.

Go to the **User interface** and under the **Toolbox** click on the link **Create MS Excel table**

Instead of typing your specimen codes in the text area below, you could select a Taxonset (provided that it has been set before (:ref:`create_taxonset`).
This table will be ready to attach to a manuscript for publication.

You can also change the way sequence information is displayed in the table by choosing between **number of bases** (displays number of bases - does not count questionmarks **?**), **accession numbers** (displays stored accession numbers instead of sequence length) or **X/-** (displays **X** if sequence is present and **-** if sequence is missing.

**Display missing sequence beginnings/ends with star(*)?:** will show search for questionmarks (?) in the beginning or end of the sequences (when displayed by number of bases) and show if the sequence misses bases in those positions with an asterisk (*). Easy then to see during laboratory phase then where sequence information might be missing for your taxa.

You may also change between comma (,) and tab-delimited table mode.

.. image:: images/create_table.png
   :align: center
   :width: 819px





.. _update_voucher:

^^^^^^^^^^^^^^
Update voucher
^^^^^^^^^^^^^^

When you **click on the code** of an already existing voucher in the **administrator interface** you will be transferred to it's **voucher information page**.

Here you may make changes to all the fields - and these will be updated after pressing **Update record**.

A changed **voucher code** will automaticly change the code in the connected fields for sequences and primer informations, so as to keep them connected.

There is also a **record history** displayed for administrators that list what changes have been made to the voucher information previously, with time and the user responsible for the changes.

.. image:: images/update_voucher.png
   :align: center
   :width: 836px



------------
Update VoSeq
------------

The easiest way to update VoSeq (that does not require new install of software or database):

* `Download the new files from github <https://github.com/carlosp420/VoSeq/tags>`_.

    * **Unpack** the new files to your **webserver directory** (htdocs, webserver, etc).
    * **Rename your old "in use" VoSeq folder** something like, 'VoSeq_old' or similar (e.g. "VoSeq-1.4.4" -> "VoSeq_1.4.4_old").
    * **Give the newly downloaded VoSeq folder the same name as the old one had** (e.g. "VoSeq-1.4.4").
    * **Copy the file "conf.php"** (in main folder) from the old version (e.g. "VoSeq_1.4.4_old") to the new version (e.g. "VoSeq_1.4.4").

Also:

* If you have used and installed blast files, make sure to copy the files **Blastdb_aliastool**, **Blast**, **Makeblastdb** and **Makembindex** (.exe for all in windows), from the old version (blast/bin folder) to the new version (same folder). Then **set permissions** to read, write, and execute on the folder "blast/bin" and its content, as well as the folder "include/blast" (e.g. ``chmod 777 -R path_to/~VoSeq_folder/blast/bin``) .

    * For Mac users it may work better to use the ``sudo chown -R _www VoSeq_folder`` command instead, since files belong to user instead of root!

* **If you have voucher photos stored**, transfer them from the old one to the new one also (in ``pictures`` folder).



--------------------------
Backup your MySQL database
--------------------------

You can make backup copies of your data by using a button in the administrator interface.
You will get all your voucher info and sequences into a SQL file. If your server dies you can easily restore your database by importing one of your backups using the **Import database** button.

.. image:: images/import_export_db.png
   :align: center
   :width: 558px


-------
Plugins
-------


^^^^^^^^^^^^^^^^^^
Yahoo! Maps plugin
^^^^^^^^^^^^^^^^^^

**VoSeq** is able to interact with Yahoo! Maps to create on-the-fly maps for vouchers when geographic coordinates are present in voucher pages.
After installing VoSeq, you can enable this capability by getting a **Yahoo! Maps API key** from them and writing them in your ``conf.php`` file:

#. Get an API key from http://developer.yahoo.com/maps/simple/
#. After filling in the required information you will be given a **Consumer Key** consisting of a long string of seemingly random characters that end with two dashes:

    * ``MwRGV2Jm1zbWNHbmnM9Y2Q9WVdrOVVHdj0yzlNQS0tJ9uc3VtZXJzZWNyZXQmeD1iMw--``

#. Remove the two dashes from the end and copy your key into the ``conf.php`` file as a value for the variable ``$yahoo_key``. Like the example below, including quotations and semicolon:

    * ``$yahoo_key = "MwRGV2Jm1zbWNHbmnM9Y2Q9WVdrOVVHdj0yzlNQS0tJ9uc3VtZXJzZWNyZXQmeD1iMw";``

#. Save the file and exit.

After doing this, VoSeq will be able to pull maps from Yahoo! whenever there is geographic information in your database. Note that you need to enter the geographic coordinates into VoSeq converted into decimal format, using the sign minus for the Southern and Western hemispheres.



.. _flickr_plugin:

^^^^^^^^^^^^^
Flickr plugin
^^^^^^^^^^^^^

**VoSeq** hosts all the specimen photos in `Flickr <http://www.flickr.com/>`_. If you have a free account you can host up to 200 photos. The Pro account allows you hosting unlimited number of photos for a yearly fee (25 USD).

#. You need to get an API key from Flickr.
#. Create and account in `Flickr <http://www.flickr.com/>`_ (if you don't own one already)
#. Go to http://nymphalidae.utu.fi/cpena/VoSeq/
#. Follow the instructions to get an **API key**, **Secret key** and **Token key**.
#. After submitting you will get your **Key**, **Secret** and **Token**. Write down those keys.
#. From a text editor software, edit the file ``conf.php`` by copying your keys in it.
#. For example [these are not real keys and will not work if you use them]:

    * ``$flickr_api_key = "2d7f59f9aaa2d5c0a2782d7f5d9083a6";``
    * ``$flickr_api_secret = "ef0def0f3d5f3f15f1";``
    * ``$flickr_api_token = "61607157718372495-f5524ead33b43129";``

#. Save and exit.

Thus, every picture that you upload into your VoSeq installation will be uploaded into your Flickr account.

.. note:: You can share your voucher photos with the Encyclopedia or Life. :ref:`sharing_photos_with_eol`




.. _blast-plugin:

^^^^^^^^^^^^
BLAST plugin
^^^^^^^^^^^^

VoSeq has `BLAST capabilities <http://en.wikipedia.org/wiki/BLAST>`_.

You can search for homologous sequences of your markers in GenBank. If you have a VoSeq installation in your work computer (or your server provider allows you to run the BLAST executable files), you do local BLASTs. For example, BLAST any or your sequences against all sequences of the same gene, or against all your sequences (full BLAST). Click on the "BLAST" icons in your voucher's pages:

.. image:: images/voseq01.png
   :align: center
   :width: 800px

You can also copy and paste any new sequence into VoSeq's **Blast new sequence** tool and see whether there are any similar sequence in your data (this tool is located on the sidebar on the right).

Remember that you need to download from NCBI the stand alone BLAST executable files and copy/install them in one of VoSeq's folders:

* In Mac OS X: when you install from the .DMG package, the executable files will be written in the folder: ``/usr/local/ncbi/blast/bin``. You just need to copy them to the right folder in VoSeq:

    * ``mkdir ~/Sites/VoSeq/blast/bin``
    * ``cp /usr/local/ncbi/blast/bin/*   ~/Sites/VoSeq/blast/bin/.``

* In Linux: ``/path/to/your/VoSeq/blast/bin/``
* In Windows: ``C:\Program Files\Apache Software Foundation\Apache2.2\htdocs\VoSeq\blast\bin\``
* It is important that the executable files are placed inside the folder **bin**.



^^^^^^^^^^^^^^^^^^^^
Integration with EOL
^^^^^^^^^^^^^^^^^^^^

#. VoSeq makes it easy to share your voucher photos with EOL. More information here :ref:`sharing_photos_with_eol`.
#. VoSeq makes automated calls to EOL's web services for pulling information on authors and date of description for species. VoSeq sends genus and species names and waits for a response. If EOL response is positive, the full species name will be included in voucher pages:

.. image:: images/authority_from_eol.png
   :align: center
   :width: 574px



.. _sharing_photos_with_eol:

^^^^^^^^^^^^^^^^^^^^^^^
Sharing Photos with EOL
^^^^^^^^^^^^^^^^^^^^^^^

VoSeq makes it easy to share your voucher photos with EOL. You can submit your best photos to EOL from VoSeq with just one click.

If you haven't done it already, you need to create an account in Flickr. Then log in to Flickr with your account and join the EOL group:

#. Go to http://www.flickr.com/groups/encyclopedia_of_life
#. Click **"Join This Group"**

Be aware that EOL requires that your photo is under any of the following licenses:

* Creative Commons Attribution (`CC-BY <http://www.flickr.com/creativecommons/>`_)
* Creative Commons Non-Commercial (`CC-BY-NC <http://www.flickr.com/creativecommons/>`_)
* Creative Commons Share-Alike (`CC-BY-SA <http://www.flickr.com/creativecommons/>`_)
* Creative Commons Non-Commercial Share Alike (`CC-BY-NC-SA <http://www.flickr.com/creativecommons/>`_)

In your VoSeq installation, you will see a button:

.. image:: images/share_with_eol.png
   :align: left
   :width: 158px

under your voucher photos. If you click this button, VoSeq will add a "machine tag" to the corresponding page in Flickr so that in can be harvested by EOL.

Thus you will be able to see your photo in EOL's pool of photos in their Flickr account http://www.flickr.com/groups/encyclopedia_of_life/pool/with/4096153224/

EOL harvests the photos quite frequently, so in one day or two, you will be able to find your photo in the respective page in EOL.





^^^^^^^^^^^^^^^^^^^^^
Integration with GBIF
^^^^^^^^^^^^^^^^^^^^^

""""""""""""""""""""""""""""""""""""""""""""""""""""""""
You can share your information hosted in VoSeq with GBIF
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

GBIF prefers data owners to use their `Integrated Publishing Toolkit (IPT) <http://www.gbif.org/informatics/infrastructure/publishing/#c889>`_. This means that you can install their IPT software to produce a resource in Darwin Core format that can be harvested by GBIF. In addition to the actual data in your VoSeq installation, IPT allows you to include a rich variety of metadata for GBIF.

VoSeq is able to produce a dump file containing all the data you own. Then you can import this file into a IPT installation and choose which types of data you want to publish via GBIF. Once you include all the metadata required by GBIF you have two choices in order to expose your data taken from `GBIF website <http://www.gbif.org/informatics/standards-and-tools/publishing-data/>`_:

* By setting up a dynamic server software:

    #. Acquire hardware with a permanent Internet connection (a regular PC is sufficient).
    #. Install data publishing software. GBIF recommends the Integrated Publishing Toolkit (IPT). You will need a web server such as Apache.
    #. Configure the software for your local data structure; this is the 'mapping' process. Please follow the documentation of your chosen publishing software for this process.
    #. Register your service with GBIF and sign the GBIF Data Sharing Agreement.
    #. Create an archive for your entire dataset:

        * This scenario doesn't require a permanent Internet connection. You simply need to create a Darwin Core Archive, upload it to a repository (for example an IPT repository installed by your GBIF Participant Node, an institutional FTP or web server, or a service like Dropbox or the Internet Archive). You then just need to register the public URL for the storage location of your archive with GBIF.



"""""""""""""""""""""""""""""""""""""""""""""""""
Create a dump file and use in an IPT installation
"""""""""""""""""""""""""""""""""""""""""""""""""

#. You can create a dump file with all the data in your VoSeq installation for submitting to GBIF. In VoSeq, on the sidebar on the right, click on **Share data with GBIF**. Save this file and open an installation of IPT.
#. In IPT, click on **Manage Resources** in the top menu, enter a name for your resource and click **Create**. Note: do not upload your file in this page (it will fail to recognize your tab delimited dump file).

.. image:: images/ipt01.png
   :align: center
   :width: 591px

#. You will be directed to your test resource page. This is when you upload the dump file generated in VoSeq. And then you are ready to add a rich variety of metadata to your resource and become a provider of information to GBIF.

.. image:: images/ipt02.png
   :align: center
   :width: 625px
