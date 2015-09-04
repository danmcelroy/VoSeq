-------
Plugins
-------

.. _google_maps_plugin:

^^^^^^^^^^^^^^^^^^
Google Maps plugin
^^^^^^^^^^^^^^^^^^

**VoSeq** is able to interact with Google Maps to create on-the-fly maps for
vouchers when geographic coordinates are present in voucher pages.
After installing VoSeq, you can enable this capability by getting a
**Google Maps API key** from them and writing them in your ``conf.json`` file:

1. Go to https://developers.google.com/maps/documentation/javascript/tutorial and
   get a **Google Maps JavaScript API** key for yourself.
2. Open your ``conf.json`` file in any text editor and write your API key by
   replacing the value in the line:

.. code-block:: javascript

   "GOOGLE_MAPS_API_KEY": "get_a_google_map_api_key"

3. Save the file, exit and VoSeq will be able to pull maps from Google whenever
   there is geographic information in your database.
   Note that you need to enter the geographic coordinates into VoSeq converted
   into decimal format, using the sign minus for the Southern and Western hemispheres.



.. _flickr_plugin:

^^^^^^^^^^^^^
Flickr plugin
^^^^^^^^^^^^^

VoSeq is able to host all the specimen photos in Flickr. If you have a free
account you can host up to 200 photos. The Pro account allows you hosting
unlimited number of photos for a yearly fee (25 USD).

You need to get `API keys from Flickr <https://www.flickr.com/services/api/keys/>`__
and place them in the ``config.json`` configuration file of VoSeq:

* Create and account in Flickr (if you don't own one already)
* Follow the instructions to get an API key and Secret key.
* After submitting you will get your Key and Secret. Write down those keys.
* Using a text editor software, edit the file ``config.json`` by copying your keys in it.

* For example [these are not real keys and will not work if you use them]:

.. code-block:: javascript

    "FLICKR_API_KEY": "2d7f59f9aaa2d5c0a2782d7f5d9083a6",
    "FLICKR_API_SECRET": "ef0def0f3d5f3f15f1"

* Save and exit.



.. _blast-plugin:

^^^^^^^^^^^^
BLAST plugin
^^^^^^^^^^^^
You can blast your sequences within VoSeq provided that the BLAST software from
NCBI is installed in your computer or server.
If you have a Ubuntu server you can easily install BLAST with the following command:

.. code-block:: shell

    sudo apt-get install ncbi-blast+

After this the BLAST tools in VoSeq should work right away.



.. _gbif:

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


