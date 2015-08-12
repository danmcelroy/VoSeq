.. _faq:

Frequently Asked Questions
==========================

Version 1.7.4 and earlier
-------------------------

Installation FAQ:
^^^^^^^^^^^^^^^^^

**Q.-** **During installation I get the error 2002: "Can't connect to local MySQL server through socket ....bla bla bla..."**

    * VoSeq is trying to connect to MySQL using a file called **socket**. This error occurs when PHP tells VoSeq to find the socket in a folder where it is not placed. This can be fixed by telling MySQL to put the socket as the file ``/tmp/mysql.sock`` and by telling PHP to find it there and not to look for it in any other folder.
    * From the installation folder of **PHP**, save the file ``php.ini-development`` in the folder ``/usr/local/lib/`` and name it ``php.ini``
    * Edit your file ``php.ini`` and look for the command ``mysql.default_socket`` and make sure it says:

        * ``mysql.default_socket = /tmp/mysql.sock``

    * Edit your MySQL installation file ``/usr/local/mysql/support-files/my-large.cnf``:

        * File parameters: modify the lines ``socket  = /var/mysqld/mysqld.sock`` to ``socket = /tmp/mysql.sock``
        * Save the file as ``/etc/my.cnf``  and ``/etc/mysql/my.cnf``

    * Restart the server and resume the installation of VoSeq.


POST request error
^^^^^^^^^^^^^^^^^^
**Q.-** ...my computer complains that "The requested resource /VoSeq_XXX/somefile.php does not allow request data with POST requests, or the amount of data provided in the request exceeds the capacity limit."?

**A.-** Open the PHP config file (see below "How to find PHP.ini") and increase the value for ``POST_max_size``, save file and restart webserver.

Execution time error
^^^^^^^^^^^^^^^^^^^^
**Q.-** ...my computer stops a VoSeq page from running due to execution timeout?

**A.-** Open the PHP config file (see below "How to find PHP.ini") and increase the value for ``max_execution_time``, save file and restart webserver.

Too many variables problem
^^^^^^^^^^^^^^^^^^^^^^^^^^
**Q.-** ...my huge taxonsets or other lists doesnt include all the values I had marked and added for them?

**A.-** PHP may have set a too low value to ``max_input_vars``. Open the PHP config file (see below "How to find PHP.ini") and increase the value for ``max_input_vars``, save file and restart webserver.

Mac permission problem
^^^^^^^^^^^^^^^^^^^^^^
**Q.-** ...if for example BLAST, storing pictures etc dont work on Mac!

**A.-** It happens specially when upgrading VoSeq. When you download a fresh copy of VoSeq form Github and copy the contents on your installation of VoSeq, it happens that all the files and folders have you as **owner**. So, VoSeq (and the Apache server) cannot write into the folders. To fix this it is necessary to set the Apache server as the **owner** of files and folder. In my MacBook the id for the Apache "user" is ``_www``. So we need to do the following to transfer ownership of files and folders of VoSeq to the server:
Open a terminal or console and use the command: ``sudo chown -R _www VoSeq_folder``. This should give the permissions to VoSeq (actually the Apache server) to do this things!

How to find php.ini and see your PHP settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Place a file named **info.php** containing ``<?php phpinfo(); ?>`` in your web server folder where you have your VoSeq folder.
Open your browser and go to that file/page (ie. **http://localhost/info.php** for win/linux or **http://127.0.0.1/~yourprivatefoldername/info.php** for mac).
This should get you the PHP config output, where you can find "Configuration File(php.ini) Path" and "Loaded Configuration File". These fields should tell you where your config file (php.ini) is located.
If these says "(none)" see below.

    * Windows - In windows the PHP configuration file (php.ini) should be found under ``C:Windows/``. If it's not there then copy the php.ini-??? to ``C:\WINDOWS`` and rename it php.ini. (??? can be dist, production or development).
    * Mac - on mac the the PHP configuration file (php.ini) should be found under ``/private/etc/`` . If no php.ini is found there but a php.ini.default is, run ``sudo cp /private/etc/php.ini.default /private/etc/php.ini`` in terminal create a php.ini file. Then restart server.
    * Linux - Open a terminal or console and type ``locate php.ini``. In my server I got this location: ``/usr/local/lib/php.ini``
