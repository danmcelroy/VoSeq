<?php
// #################################################################################
// #################################################################################
// Voseq admin/processfile.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: process picture file uploaded by user
//
// #################################################################################


// #################################################################################
// Section: include functions
// #################################################################################

//check admin login session
include'../login/auth-admin.php';
// includes
ob_start();//Hook output buffer - disallows web printing of file info...
include'../conf.php';
ob_end_clean();//Clear output buffer//includes
#include'../login/redirect.html';

include 'admarkup-functions.php';
include '../includes/smart_resize_image.php';

if( $photos_repository == "flickr" ) {
	require_once'../api/phpFlickr/phpFlickr.php';

	// create an api
	$f = new phpFlickr($flickr_api_key, $flickr_api_secret);
	$f->setToken($flickr_api_token);
}

// array to keep track of the photo ids as they're uploaded
$photo_ids = array();


function strip_ext($name) {
	$ext = strrchr($name, '.');
	if($ext !== false) {
		$name = substr($name, 0, -strlen($ext));
	}
   return $name;
}


# Support for several photos per voucher
# photo details will be appended and divided using "palote" |
# (p_, thumbnail, flickr_id, voucherImage, code)
function update_photo_details($p_, $my_url, $photo_id, $my_voucherImage, $code) { 
	# See if we are already have a photo for this voucher
	$query  = "SELECT id, thumbnail, flickr_id, voucherImage FROM ";
	$query .= $p_ . "vouchers WHERE code='" . $code . "'";

	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	# if records present append to MySQL
	if (mysql_num_rows($result) > 0) {
		while ($row = mysql_fetch_object($result)) {
			$photos = $row;
			
			# remove na.gif if necessary
			if( $photos->thumbnail == "na.gif" && $photos->voucherImage == "na.gif" ) {
				$photos->thumbnail = "";
				$photos->voucherImage = "";
			}
			$photos->thumbnail .= "|" . $my_url;
			if( trim($photo_id) == "" ) {
				$photo_id = "";
			}
			$photos->flickr_id .= "|" . $photo_id;
			$photos->voucherImage .= "|" . $my_voucherImage;
		}
		$query  = "UPDATE ". $p_ . "vouchers set timestamp=now(), ";
		$query .= "thumbnail=\"$photos->thumbnail\", ";
		$query .= "flickr_id=\"$photos->flickr_id\", ";
		$query .= "voucherImage=\"$photos->voucherImage\" ";
		$query .= "where code=\"$code\""; 
		mysql_query($query) or die("Error in query: $query. " . mysql_error());
	}
	# no records present, fresh upload of data to MySQL
	else {
		$query  = "UPDATE ". $p_ . "vouchers set timestamp=now(), ";
		$query .= "thumbnail=\"$my_url\", flickr_id=\"$photo_id\", ";
		$query .= "voucherImage=\"$my_voucherImage\" where code=\"$code\""; 
		mysql_query($query) or die("Error in query: $query. " . mysql_error());
	}
}


// to indicate this is an administrator page
$admin = true;

// check for record ID
if ((!isset($_GET['code']) || trim($_GET['code']) == '')) {
	die('Missing record ID!');
	}
else {
	$code = $_GET['code'];
 	$a = "{$_FILES['userfile']['name']}";
	$b = trim($a);
	if ($b == '') {
		// process title
		$title = "$config_sitename -Upload failed";

		// print html headers
		include_once('../includes/header.php');

		// print navegation bar
		admin_nav();

		// begin HTML page content
		echo "<div id=\"content_narrow\">";
		echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
				<tr><td valign=\"top\">";
		echo "<img src=\"../images/warning.png\" alt=\"\"> ";
		echo "File <b>did not</b> upload.<br />Check the file size. ";
		echo "File must be less than 4MB. Or maybe your filename is not correct.<br />";
		echo "Another posibility is that your <code>php.ini</code> file has been set ";
		echo "to accept only small files and you might need to edit it.";
		echo "</td>";
		echo "<td class=\"sidebar\" valign=\"top\">";
		admin_make_sidebar();
		echo "</td>";
		echo "</tr>
				</table> <!-- end super table -->
				</div> <!-- end content -->";
		make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
		echo "\n</body>\n</html>";
	}
	else {
		// open database connections
		@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect');
		mysql_select_db($db) or die ('Unable to select database');
		if( function_exists(mysql_set_charset) ) {
			mysql_set_charset("utf8");
		}

		// generate and execute query
		$query = "SELECT genus, species, subspecies, family, subfamily, tribe, subtribe, country, specificLocality, publishedIn, notes, voucherImage, latitude, longitude FROM ". $p_ . "vouchers WHERE code=\"$code\"";
		$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
			
		// if records present
		if (mysql_num_rows($result) > 0) {
			while ($row = mysql_fetch_object($result)) {
				$genus = $row->genus;
				$species = $row->species;
				$subspecies = $row->subspecies;
				$family = $row->family;
				$subfamily = $row->subfamily;
				$tribe = $row->tribe;
				$subtribe = $row->subtribe;
				$latitude = $row->latitude;
				$longitude = $row->longitude;
				if( $row->country != "" ) {
					$country = "$row->country. ";
				}
				else {
					$country = "";
				}

				if( $row->specificLocality != "" ) {
					$specificLocality = "$row->specificLocality, ";
				}
				else {
					$specificLocality = "";
				}

				if( $row->publishedIn != "" ) {
					$publishedIn = "$row->publishedIn, ";
				}
				else {
					$publishedIn = "";
				}

				if( $row->notes != "" ) {
					$notes = "$row->notes, ";
				}
				else {
					$notes = "";
				}
			}
		}
		
		// process title
		$title = "$config_sitename - Picture uploaded ";

		// print html headers
		include_once('../includes/header.php');

		// print navegation bar
		admin_nav();

		// begin HTML page content
		echo "<div id=\"content_narrow\">";
		echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
				<tr><td valign=\"top\">";
		
		$item = $_FILES['userfile']['tmp_name'];
		$extension = substr($item, - strlen(PHOTO_EXTENSION));
		print "Uploading $item...\n";

		if( $photos_repository == "flickr" ) {
			$photo_id = $f->sync_upload($item, "$code $genus $species $subspecies", "$country $specificLocality $publishedIn $notes <a href=\"$base_url/story.php?code=$code\">see in our database</a>", "$country,$family,$subfamily,$tribe,$subtribe,$genus,$species,$subspecies");
			$info = $f->photos_getInfo($photo_id);
			$sizes = $f->photos_getSizes($photo_id);
			$my_voucherImage = $info['photo']['urls']['url'][0]['_content'];
			$status = $f->photos_geo_setLocation($photo_id, $latitude, $longitude, "3");
		
			/*** create thumbnails ***/
			foreach( $sizes as $item) {
				foreach($item as $k => $v) {
					if($k == "width" && $v == "240") {
						$my_url = $item['source'];
					}
					elseif($k == 'label' && $v == 'Small') {
						$my_url = $item['source'];
					}
					elseif($k == 'label' && $v == 'Thumbnail') {
						$my_url = $item['source'];
					}
				}
			}

			# Support for several photos per voucher
			# photo details will be appended and divided using "palote" |
			# (thumbnail, flickr_id, voucherImage, code)
			update_photo_details($p_, $my_url, $photo_id, $my_voucherImage, $code); 
		}
		else {  // $photos_repository == "local" 
			$dest_folder = $local_folder . "/pictures";

			$filename = $_FILES['userfile']['name'];
			// check that file doesnot exist to avoid overwritting
			$i = 1;
			while ( file_exists("$dest_folder/$filename") ) {
				$filename = $_FILES['userfile']['name'];
				$filename = $i . "_" . $filename;
				$i = $i + 1;
			}
			echo " as file name: <b>$filename</b>";
			$tmp_name = $_FILES['userfile']['tmp_name'];
			move_uploaded_file($tmp_name, "$dest_folder/$filename");

			if( !function_exists('imagecreatefromjpeg') ) {
				echo "<p> <img src='../images/warning.png' />Warning: You need to ";
				echo "install the library GD in your computer in ";
				echo "order to create thumbnails for images.</p>";
				echo "<p>If you compiled PHP from source, compile it with support for GD.";
				echo "<br />I will proceed but it would be best if you install GD! ";
				echo "See more info <a href='http://www.php.net/manual/en/image.installation.php'>";
				echo "here</a>.</p>";
				$thumbnail_name = $base_url . "/pictures/" . $filename;
			}	
			else {
				$thumbnail_name = "$dest_folder/thumbnails/tn_$filename";
				smart_resize_image("$dest_folder/$filename",
								200,
								0,
								true,
								$thumbnail_name,
								$delete_original=false,
								false
								);
				$thumbnail_name = $base_url . "/pictures/thumbnails/tn_" . $filename;
			}
			$filename = $base_url . "/pictures/" . $filename;

			# Support for several photos per voucher
			# photo details will be appended and divided using "palote" |
			# (thumbnail, flickr_id, voucherImage, code)
			### $query = "UPDATE ". $p_ . "vouchers set timestamp=now(), thumbnail=\"$thumbnail_name\", flickr_id=null, voucherImage=\"$filename\" where code=\"$code\""; 
			update_photo_details($p_, $thumbnail_name, NULL, $filename, $code); 
		}



		echo "\n";
			
		echo "<br /><img src=\"images/success.png\" alt=\"\" /> File is valid and was sucessfuly uploaded. " . $subject . " (" . $_FILES['userfile']['size'] . ")";
	
		
		// show links to add sequences for this same record
		?>
		<br /><br />
		Do you want to:
		<ol>
			<li>Enter sequences for record of code <b><?php echo "$code"; ?></b>:
				<?php 
				if( $mask_url == "true" ) {
					echo "<a href='" .$base_url . "/home.php' onclick=\"return redirect('listseq.php?code=$code')\">Add Sequences</a></li>";
					echo "<li><a href='" .$base_url . "/home.php' onclick=\"return redirect('admin.php')\">Go back to the main menu</a>.</li>";
				}
				else {
					echo "<a href='" .$base_url . "/admin/listseq.php?code=$code'>Add Sequences</a></li>";
					echo "<li><a href='" .$base_url . "/admin/admin.php'>Go back to the main menu</a>.</li>";
				}
		echo "</ol>";
		echo "</td>";
		echo "<td class=\"sidebar\" valign=\"top\">";
		admin_make_sidebar();
		echo "</td>";
		echo "</tr>
				</table> <!-- end super table -->
				</div> <!-- end content -->";
		make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
		echo "\n</body>\n</html>";
			}
		}	
	
?>
