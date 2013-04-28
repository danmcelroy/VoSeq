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

	// table definition for voucher images changed to accomodate several strings separated by "|".
	// version 1.5.0
	$query = "SELECT * FROM information_schema.COLUMNS WHERE 
				TABLE_SCHEMA = '" . $db . "' 
				AND TABLE_NAME = '" . $p_ . "vouchers' 
				AND COLUMN_NAME = 'flickr_id'";

	$query  = "SELECT data_type FROM information_schema.COLUMNS WHERE
				TABLE_SCHEMA = '" . $db . "' 
				AND TABLE_NAME = '" . $p_ . "vouchers'
				AND COLUMN_NAME = 'thumbnail'";
	$result = mysql_query($query) or die ("Error in query: $query. " . mysql_error());
	if( mysql_num_rows($result) > 0 ) {
		while( $row = mysql_fetch_object($result) ) {
			if( $row->data_type == "varchar" ) {
				$query = "ALTER TABLE ". $p_ . "vouchers MODIFY voucherImage text default null";
				mysql_query($query) or die ("Error in query: $query. " . mysql_error());
				$query = "ALTER TABLE ". $p_ . "vouchers MODIFY flickr_id text default null";
				mysql_query($query) or die ("Error in query: $query. " . mysql_error());
				$query = "ALTER TABLE ". $p_ . "vouchers MODIFY thumbnail text default null";
				mysql_query($query) or die ("Error in query: $query. " . mysql_error());
			}
		}
	}
	
	
	// genes -> genetype, prot_code, intron, aligned, genetic_code - version 1.4.5
	$new_fields = array("genetype", "prot_code", "intron", "aligned", "genetic_code");
	foreach ($new_fields as $new_field) {
		$query = "SELECT * FROM information_schema.COLUMNS WHERE 
					TABLE_SCHEMA = '" . $db . "' 
					AND TABLE_NAME = '" . $p_ . "genes' 
					AND COLUMN_NAME = '$new_field'";
		$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
		if( mysql_num_rows($result) < 1 ) {
			# Upgrading table sequences
			$query = "alter table ". $p_ . "genes add column $new_field";
			if ($new_field == 'genetype'){$query .= " varchar(255) default null";}
			if ($new_field == 'prot_code' || $new_field == 'aligned'){$query .= " ENUM('yes','no','notset') default 'notset'";}
			if ($new_field == 'intron'){$query .= " varchar(255) default null";}
			if ($new_field == 'genetic_code'){$query .= " INT(3) default null";}
			$query .= " after timestamp";
			mysql_query($query) or die ("Error in query: $query. " . mysql_error());
		}
	}
	
	// genesets - from v.1.7.0
	$query = "SELECT *
		FROM information_schema.tables 
		WHERE table_schema = '" . $db . "' 
		AND table_name = '" . $p_ . "genesets'";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	if( mysql_num_rows($result) < 1 ) {
		$query = "CREATE TABLE `" . $db . "` . `" . $p_ . "genesets` (
			`geneset_name` varchar(75) DEFAULT NULL,
			`geneset_creator` varchar(75) DEFAULT NULL,
			`geneset_description` varchar(100) DEFAULT NULL,
			`geneset_list` text,
			`geneset_id` int(11) NOT NULL AUTO_INCREMENT,
			PRIMARY KEY (`geneset_id`),
			UNIQUE KEY `geneset_id_UNIQUE` (`geneset_id`)
			) ENGINE=MyISAM AUTO_INCREMENT=6 DEFAULT CHARSET=latin1;";
		mysql_query($query) or die("Error in query: $query. " . mysql_error());
	}
}
	
?>
