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
