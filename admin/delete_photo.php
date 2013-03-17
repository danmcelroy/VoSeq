<?php



# Delete photo. It needs voucherImage and thumbnail names
$voucherImage = $_POST['voucherImage'];
$thumbnail = $_POST['thumbnail'];
preg_match("/http:\/\/.+(\/.+\/)$/", $voucherImage, $match);

if( count($match) > 0 ) {
	$flickr_id = str_replace("/", "", $match[1]);
}
else {
	$flickr_id = " ";
}

if( $voucherImage != "" && $thumbnail != "" ) {
	ob_start();
	include('../conf.php');
	ob_end_clean();

	// open database connection
	$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
	//select database
	mysql_select_db($db) or die ('Unable to content');
	if( function_exists(mysql_set_charset) ) {
		mysql_set_charset("utf8");
	}

	$query  = "SELECT id, voucherImage, flickr_id, thumbnail FROM ". $p_;
	$query .= "vouchers WHERE voucherImage ";
	$query .= "LIKE '%" . $voucherImage . "%'";
	$result = mysql_query($query) or die("Error in query: $query ". mysql_error());
	if( mysql_num_rows($result) > 0 ) {
		while( $row = mysql_fetch_object($result) ) {
			$voucherImage = "|" . $voucherImage;
			$v = str_replace($voucherImage, "", $row->voucherImage);

			$flickr_id = "|" . $flickr_id;
			$f = str_replace($flickr_id, "", $row->flickr_id);

			$thumbnail = "|" . $thumbnail;
			$t = str_replace($thumbnail, "", $row->thumbnail);

			$q  = "UPDATE " . $p_ . "vouchers SET ";
			$q .= "voucherImage='$v', flickr_id='$f', thumbnail='$t' ";
			$q .= "WHERE id=$row->id";	
			mysql_query($q) or die("Error in query: $q ". mysql_error());
		}
	}
}
?>
