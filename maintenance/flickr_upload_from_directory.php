#!/usr/local/bin/php
<?php
error_reporting(0);

if (count($argv) < 2) {
    echo "\nError, please enter the folder containing the pictures to upload";
    echo " as argument.\n\n";
    echo "\t" . "php flickr_upload_from_directory.php /home/user/Desktop/photos\n";
    echo "\nExiting. Nothing done.\n";
    exit(0);
}


# check if conf file is present
if (file_exists('../conf.php')) {
    include_once('../conf.php');
}
else {
    echo "\nError, couldn't find the configuration file ``conf.php`` in";
    echo " VoSeq's main folder.";
    echo "\nExiting. Nothing done.\n";
    exit(0);
}


# check for Flickr API Keys
$error = "";
if ($photos_repository == "flickr") {
    echo "\nI will upload all the fotos to Flickr, as stated in your";
    echo " ``conf.php`` file.\n";
    if (!isset($flickr_api_key)) {
        $error = "true";
    }
    elseif (isset($flickr_api_key) and strlen($flickr_api_key) < 1) {
        $error = "true";
    }
    elseif (!isset($flickr_api_secret)) {
        $error = "true";
    }
    elseif (isset($flickr_api_secret) and strlen($flickr_api_secret) < 1) {
        $error = "true";
    }
    elseif (!isset($flickr_api_token)) {
        $error = "true";
    }
    elseif (isset($flickr_api_token) and strlen($flickr_api_token) < 1) {
        $error = "true";
    }

    if ($error == "true") {
        $error = "\nError, you need to get the Flickr API keys and include them";
        $error .= " in the ``conf.php`` file\n";
        $error .= "See here for more info: ";
        $error .= "http://nymphalidae.utu.fi/cpena/VoSeq_docu.html#flickr-plugin\n";
        $error .= "\nExiting. Nothing done.\n";
        echo $error;
        exit(0);
    }
}
elseif ($photos_repository == "local") {
    echo "\nI will upload all the fotos to the local instalation of VoSeq";
    echo " and all the files will be kept into the local hardisk,";
    echo " as stated in your ``conf.php`` file.\n";
    echo "\nNoting will be sent to Flickr.\n";
}
else {
    $error = "\nError, you need to specify in your ``conf.php`` file\n";
    $error .= " whether the files will be sent to Flickr or kept only in";
    $error .= " your harddrive.\n";
    $error .= "\nExiting. Nothing done.\n";
    echo $error;
    exit(0);
}

echo "\n===========================\n";
echo "\nTo upload voucher pics from a directory you need to make sure that the";
echo "file name is the same as the voucher code.\n";
echo "Also make sure that the voucher info is already in the database.\n";
echo "\n===========================\n";

define('UPLOAD_DIRECTORY', $argv[1]);
if (file_exists(UPLOAD_DIRECTORY)) {
    echo "\nI will upload the pictures from this folder: ";
    echo UPLOAD_DIRECTORY . "\n";
}
else {
    echo "\nError, couldn't find the folder ``". UPLOAD_DIRECTORY . "``";
    echo "\nExiting. Nothing done.\n";
    exit(0);
}

require_once('../api/phpFlickr/phpFlickr.php');

# define('PHOTO_EXTENSION', '.png');

# create api
$f = new phpFlickr($flickr_api_key, $flickr_api_secret);
$f->setToken($flickr_api_token);

# create a DirectoryIterator (part of the Standard PHP Library)
$di = new DirectoryIterator(UPLOAD_DIRECTORY);

# open database connections
$mysqli = mysqli_init();
if (!$mysqli) {
    die('mysqli_init failed');
}
if (!$mysqli->options(MYSQLI_INIT_COMMAND, 'SET AUTOCOMMIT = 0')) {
    die('Setting MYSQLI_INIT_COMMAND failed');
}

if (!$mysqli->options(MYSQLI_OPT_CONNECT_TIMEOUT, 5)) {
    die('Setting MYSQLI_OPT_CONNECT_TIMEOUT failed');
}

if (!$mysqli->real_connect($host, $user, $pass, $db)) {
    die('Unable to connect to MySQL ('. mysqli_connect_errno() . ') '
            . mysqli_connect_error());
}

echo "Success. Connected to MySQL... " . $mysqli->host_info . "\n";
$mysqli->query("set names utf8") or die("Error in query: " . $mysqli->error());

$photos = array();
foreach($di as $file) {
    preg_match("/^(.+)\.png\$/", $file, $matches);
    if( $matches ) {
        $photos[] = $matches[1];
    }
}

$mysqli->close();
// upload pictures in the directory
foreach( $photos as $item ) {
    $query = "SELECT id, code, genus, species, subspecies, family, subfamily, tribe, subtribe, country, specificLocality, publishedIn, notes, voucherImage, latitude, longitude FROM ". $p_ . "vouchers WHERE code = \"$item\"";
    $result = mysql_query($query) or die("Error in query: $query. " . mysql_error());

    while( $row = mysql_fetch_object($result) ) {
        $code = $row->code;
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
        
    $file = UPLOAD_DIRECTORY . "/" . $item . PHOTO_EXTENSION;

    $photo_id = $f->sync_upload($file, "$code $genus $species $subspecies", "$country $specificLocality $publishedIn $notes", "$country,$family,$subfamily,$tribe,$subtribe,$genus,$species,$subspecies");

    $info = $f->photos_getInfo($photo_id);
    $my_voucherImage = $info['photo']['urls']['url'][0]['_content'];
    $status = $f->photos_geo_setLocation($photo_id, $latitude, $longitude, "3");
    $sizes = $f->photos_getSizes($photo_id);
    
    /*** create thumbnails ***/
    foreach( $sizes as $i) {
        foreach($i as $k => $v) {
            if($k == "label" && $v == "Small") {
                $my_url = $i['source'];
            }
        }
    }
    $query = "UPDATE ". $p_ . "vouchers set timestamp=now(), thumbnail=\"$my_url\", flickr_id=\"$photo_id\", voucherImage=\"$my_voucherImage\" where code=\"$item\""; 
    echo $query . ";\n";
#mysql_query($query) or die("Error in query: $query. " . mysql_error());
}


?>
