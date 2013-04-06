<?php
include("../conf.php");

@$connection = mysql_connect($host, $user, $pass) or die("Unable to connect to MySQL");
mysql_select_db($db) or die ("Unable to select database");
mysql_query("set names utf8") or die ("Error in query: $query. ". mysql_error());


$geneCode = mysql_real_escape_string($_POST['id']);
$query = "DELETE FROM ". $p_ . "genes WHERE geneCode='$geneCode'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());

$query = "DELETE FROM ". $p_ . "primers WHERE geneCode='$geneCode'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());

$query = "DELETE FROM ". $p_ . "sequences WHERE geneCode='$geneCode'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());


$query = "SELECT geneCode FROM ". $p_ . "genes WHERE geneCode = '$geneCode'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());
if( mysql_num_rows($result) > 0 ) {
	echo "fail";
	exit(0);
}
else {
	echo "ok";
	exit(0);
}


?>
