<?php
include("../conf.php");

@$connection = mysql_connect($host, $user, $pass) or die("Unable to connect to MySQL");
mysql_select_db($db) or die ("Unable to select database");
mysql_query("set names utf8") or die ("Error in query: $query. ". mysql_error());


$id = mysql_real_escape_string($_POST['id']);

$query = "SELECT code, geneCode FROM ". $p_ . "sequences WHERE id='$id'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());
$row    = mysql_fetch_object($result);

$code = $row->code;
$geneCode = $row->geneCode;

$query = "DELETE FROM ". $p_ . "primers WHERE code='$code' AND geneCode='$geneCode'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());

$query = "DELETE FROM ". $p_ . "sequences WHERE code='$code' AND geneCode='$geneCode'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());


$query = "SELECT id FROM ". $p_ . "sequences WHERE code='$code' AND geneCode='$geneCode'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());
if( mysql_num_rows($result) > 0 ) {
	echo "fail"; //_id=$id code=$code geneCode=$geneCode";
	exit(0);
}
else {
	echo "ok";
	exit(0);
}


?>
