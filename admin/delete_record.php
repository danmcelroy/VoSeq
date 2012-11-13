<?php
include("conf.php");

@$connection = mysql_connect($host, $user, $pass) or die("Unable to connect to MySQL");
mysql_select_db($db) or die ("Unable to select database");
mysql_query("set names utf8") or die ("Error in query: $query. ". mysql_error());

$id = mysql_real_escape_string($_POST['id']);

$query = "DELETE FROM ". $p_ . "refs WHERE id='$id'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());

$query = "SELECT id FROM ". $p_ . "refs WHERE id='$id'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());

if( mysql_num_rows($result) > 0 ) {
	echo "fail";
}
else {
	echo "ok";
}

?>
