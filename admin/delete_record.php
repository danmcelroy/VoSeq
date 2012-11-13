<?php
include("../conf.php");

@$connection = mysql_connect($host, $user, $pass) or die("Unable to connect to MySQL");
mysql_select_db($db) or die ("Unable to select database");
mysql_query("set names utf8") or die ("Error in query: $query. ". mysql_error());

$code = mysql_real_escape_string($_POST['id']);

$query = "DELETE FROM ". $p_ . "vouchers WHERE code='$code'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());

$query = "DELETE FROM ". $p_ . "primers WHERE code='$code'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());

$query = "DELETE FROM ". $p_ . "sequences WHERE code='$code'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());

$query = "SELECT taxonset_id, taxonset_list FROM ". $p_ . "taxonsets WHERE taxonset_list like '%$code%'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());
if( mysql_num_rows($result) > 0 ) {
	while( $row = mysql_fetch_object($result) ) {
		$taxon_list_array = explode(",", $row->taxonset_list);
		$new_taxon_list = array();
		foreach($taxon_list_array as $item) {
			if( strtoupper($item) != strtoupper($code) ) {
				$new_taxon_list[] = $item;
			}
		}
		$new_taxon_list = implode(",", $new_taxon_list);
		$query2 = "update ". $p_ . "taxonsets set taxonset_list='$new_taxon_list' where taxonset_id='$row->taxonset_id'";
		$result2 = mysql_query($query2) or die("Error in query: $query2.".mysql_error());
	}
}

$query = "SELECT taxonset_id, taxonset_list FROM ". $p_ . "taxonsets WHERE taxonset_list like '%$code%'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());
if( mysql_num_rows($result) > 0 ) {
	echo "fail";
	exit(0);
}
else {
	$query = "SELECT code FROM ". $p_ . "vouchers WHERE code = '$code'";
	$result = mysql_query($query) or die("Error in query: $query.".mysql_error());
	if( mysql_num_rows($result) > 0 ) {
		echo "fail";
		exit(0);
	}
	else {
		echo "ok";
		exit(0);
	}
}

?>
