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

$query = "SELECT geneset_id, geneset_list FROM ". $p_ . "genesets WHERE geneset_list like '%$geneCode%'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());
if( mysql_num_rows($result) > 0 ) {
	while( $row = mysql_fetch_object($result) ) {
		$gene_list_array = explode(",", $row->geneset_list);
		$new_gene_list = array();
		foreach($gene_list_array as $item) {
			if( strtoupper($item) != strtoupper($geneCode) ) {
				$new_gene_list[] = $item;
			}
		}
		$new_gene_list = implode(",", $new_gene_list);
		$query2 = "update ". $p_ . "genesets set geneset_list='$new_gene_list' where geneset_id='$row->geneset_id'";
		$result2 = mysql_query($query2) or die("Error in query: $query2.".mysql_error());
	}
}


$query = "SELECT geneset_id, geneset_list FROM ". $p_ . "genesets WHERE geneset_list like '%$geneCode%'";
$result = mysql_query($query) or die("Error in query: $query.".mysql_error());
if( mysql_num_rows($result) > 0 ) {
	echo "fail";
	exit(0);
}
else {
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
}

?>
