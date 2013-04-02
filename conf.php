<?php
error_reporting (E_ALL ^ E_NOTICE);

session_start();



$intro_msg = "<h2>This is <b>VoSeq</b>. VoSeq is a database to store voucher and sequence data.
Please send all bug complaints to <br />Carlos Peña (<i>mycalesis@gmail.com</i>) or<br />
Tobias Malm (<i>tobias.malm@uef.fi</i>)&nbsp;&nbsp;</h2>";


# prefix for tables
$p_ = '';



$mask_url = 'false';

$host = 'localhost';
$user = 'root';
$pass = 'mysqlboriska';
$db   = 'filemaker';

$config_sitename = 'mybutts';
$version = "1.3.6";
$currentTemplate = 'templates/mytemplate';

$base_url = 'http://localhost/VoSeq';
$local_folder = '/home/carlosp420/data/VoSeq';

$date_timezone = "Europe/Helsinki"; // php5
$php_version = "5";

$flickr_api_key = "40f181d15c560e9082b6f342672cdd85";
$flickr_api_secret = "1522b3e873a919d3";
$flickr_api_token = "72157629518198099-483c963eb77e6767";

$yahoo_key = "dj0yJmk9OWg2VmNPaHhvNWZFJmQ9WVdrOVRtNXZabkZoTnpZbWNHbzlNVGN5TkRBeE5qWTJNZy0tJnM9Y29uc3VtZXJzZWNyZXQmeD01YQ"; # You need to get an API key for Yahoo! Maps, it's free! This is a test API key


// What repository you want to use to store voucher photos?
$photos_repository = 'local';

?>
