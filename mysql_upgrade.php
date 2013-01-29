<?php

// #################################################################################
// Section: Upgrade MySQL tables if needed
// #################################################################################

function mysql_upgrade($db, $p_) {
	// vouchers -> determinedBy
	$query = "SELECT * FROM information_schema.COLUMNS WHERE 
				TABLE_SCHEMA = '" . $db . "' 
				AND TABLE_NAME = '" . $p_ . "vouchers' 
				AND COLUMN_NAME = 'determinedBy'";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	if( mysql_num_rows($result) < 1 ) {
		# Upgrading table vouchers
		$query = "alter table ". $p_ . "vouchers add column determinedBy varchar(255)";
		$query .= " default null after timestamp";
		mysql_query($query) or die ("Error in query: $query. " . mysql_error());
	}
	

	// vouchers -> auctor
	$query = "SELECT * FROM information_schema.COLUMNS WHERE 
				TABLE_SCHEMA = '" . $db . "' 
				AND TABLE_NAME = '" . $p_ . "vouchers' 
				AND COLUMN_NAME = 'auctor'";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	if( mysql_num_rows($result) < 1 ) {
		# Upgrading table vouchers
		$query = "alter table ". $p_ . "vouchers add column auctor varchar(255)";
		$query .= " default null after timestamp";
		mysql_query($query) or die ("Error in query: $query. " . mysql_error());
	}
	

	// vouchers -> notes	
	$query = "SELECT * FROM information_schema.COLUMNS WHERE 
				TABLE_SCHEMA = '" . $db . "' 
				AND TABLE_NAME = '" . $p_ . "vouchers' 
				AND COLUMN_NAME = 'notes'";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	if( mysql_num_rows($result) < 1 ) {
		# Upgrading table sequences
		$query = "alter table ". $p_ . "vouchers add column notes varchar(255)";
		$query .= " default null after timestamp";
		mysql_query($query) or die ("Error in query: $query. " . mysql_error());
	}

	
	// vouchers -> flickr_id	
	$query = "SELECT * FROM information_schema.COLUMNS WHERE 
				TABLE_SCHEMA = '" . $db . "' 
				AND TABLE_NAME = '" . $p_ . "vouchers' 
				AND COLUMN_NAME = 'flickr_id'";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	if( mysql_num_rows($result) < 1 ) {
		# Upgrading table vouchers
		$query = "alter table ". $p_ . "vouchers add column flickr_id varchar(255)";
		$query .= " default null after timestamp";
		mysql_query($query) or die ("Error in query: $query. " . mysql_error());
	}

	
	// sequences -> notes	
	$query = "SELECT * FROM information_schema.COLUMNS WHERE 
				TABLE_SCHEMA = '" . $db . "' 
				AND TABLE_NAME = '" . $p_ . "sequences' 
				AND COLUMN_NAME = 'notes'";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	if( mysql_num_rows($result) < 1 ) {
		# Upgrading table sequences
		$query = "alter table ". $p_ . "sequences add column notes varchar(255)";
		$query .= " default null after timestamp";
		mysql_query($query) or die ("Error in query: $query. " . mysql_error());
	}
}
	
?>
