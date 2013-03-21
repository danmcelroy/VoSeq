<?php
// #################################################################################
// #################################################################################
// Voseq admin/process_upload_sequences.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: proccess sequences to upload
//
// #################################################################################


// #################################################################################
// Section: include functions
// #################################################################################

//check admin login session
include'../login/auth-admin.php';
// includes
include 'admarkup-functions.php';
ob_start();//Hook output buffer - disallows web printing of file info...
include'../conf.php';
ob_end_clean();//Clear output buffer//includes
include '../functions.php';
include '../includes/validate_coords.php';

set_time_limit(120);

error_reporting (E_ALL ^ E_NOTICE);

// to indicate this is an administrator page
$admin = true;

$charset_count = array();
unset($errorList, $field_array, $code_list, $code_array, $geneCode_array);


// #################################################################################
// Section: sanitize strings
// #################################################################################
function clean_item ($item) {
	$item = stripslashes($item);
	$item = str_replace("'", "", $item);
	$item = str_replace('"', "", $item);
	$item = str_replace(',', "", $item);
	$item = preg_replace('/^\s+/', '', $item);
	$item = preg_replace('/\s+$/', '', $item);
	$item = strtolower($item);
	return $item;
}


// #################################################################################
// Section: to show errors
// #################################################################################
function show_errors($se_in) {
	// error found
	// print navegation bar
	admin_nav();
	// begin HTML page content
	echo "<div id=\"content_narrow\">";
	echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
			<tr><td valign=\"top\">";
	// print as list
	echo "<img src=\"../images/warning.png\" alt=\"\"> The following errors were encountered:";
	echo '<br>';
	echo '<ul>';
		$se_in = array_unique($se_in);
		$se_in[] = "</br>Nothing added to db, Please revise your data!"; 
	foreach($se_in AS $item) {
		echo "$item</br>";
	}
	echo "</td>";
	
	echo "<td class=\"sidebar\" valign=\"top\">";
	admin_make_sidebar();
	echo "</td>";
	echo "</tr>
			</table> <!-- end super table -->
			</div> <!-- end content -->";
	//make footer
	make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
	?></body></html><?php
}


$input_type = $_POST['seqorvouch'];

//fix field-value comparison list
if ($input_type == 'vouch') {
	$field_values = array("code","order","family","subfamily","tribe","subtribe","genus","species","subspecies","auctor","hostorg","typespecies","country","locality","longitude","latitude","altitude","collector","coll.date","vouchercode","voucher","voucherlocality","determined.by","sex","extraction","extractiontube","extractor","extr.date","publ.in", "notes");
}
else {
	$field_values = array("code","genecode","sequences","laborator","accession","primer1","primer2","primer3","primer4","primer5","primer6","notes");
}


$field_value_count = '0'; //for counting how many fields we have


//fix input into array with lines representing vouchers -> first line = field names
$raw_voucher_upload = str_replace("\r\n", "\n", $_POST['input_data']);
$raw_voucher_upload = str_replace("\n\n", "\n", $raw_voucher_upload);
//$raw_voucher_upload = str_replace("	", "','", $raw_voucher_upload);
$raw_voucher_upload = rtrim($raw_voucher_upload, "\n\r\0 \x0B");

$lines = explode("\n", $raw_voucher_upload);



# trim white spaces from values
$tmp_line = array();
foreach($lines as $line) {
	$line = explode("	", $line);

	$tmp = array();
	foreach($line as $item) {
		$item = trim($item);
		$tmp[] = $item;
	}
	$line = implode("	", $tmp);
	unset($tmp);

	$tmp_line[] = $line;
}
$lines = $tmp_line;



//get first line of array - the field header line into its own array
$fields = explode("	", strtolower(array_shift($lines)));

///find out if we have a code field at all
if ( !in_array('code', $fields)) {
	$errorList[] = "There is no <b>Code</b> field at all
					</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;can't proceed without unique codes!" ;
}
if ( $input_type == 'seq') {
	///find out if we have a geneCode field at all
	if ( ! in_array('genecode', $fields)){
		$errorList[] = "There is no <b>Genecode</b> field at all
						</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;no use to proceed without genecodes!" ;
	}
}
// check for duplicate headers
if (count($fields) != count(array_unique)) {
	foreach ($fields as $item){
		$num_field = array_keys($fields, $item);
		if (count($num_field) > 1){
			$errorList[] = "Duplicate field headers: there is two columns named: <b>$item</b>" ;
		}
	}
	unset($item);
}
//testing if field headers are legal and adds them with the right name to a for-use array ($field_array)
if ( $input_type == 'vouch') {
	foreach($fields AS $item) {
		$field_value_count = $field_value_count + 1;
		if (in_array($item, $field_values)){
			if ($item == "order") { $field_array[] = "orden" ;}
			elseif ($item == "typespecies") { $field_array[] = "typeSpecies" ;}
			elseif ($item == "locality") { $field_array[] = "specificLocality" ;}
			elseif ($item == "coll.date") { $field_array[] = "dateCollection" ;}
			elseif ($item == "voucherlocality") { $field_array[] = "voucherLocality" ;}
			elseif ($item == "vouchercode") { $field_array[] = "voucherCode" ;}
			elseif ($item == "extr.date") { $field_array[] = "dateExtraction" ;}
			elseif ($item == "extractiontube") { $field_array[] = "extractionTube" ;}
			elseif ($item == "publ.in") { $field_array[] = "publishedIn" ;}
			elseif ($item == "determined.by") { $field_array[] = "determinedBy" ;}
			else { $field_array[] = $item ;}
			}
		else {
			$errorList[] = "Invalid field header: <b>$item</b> in column <b>$field_value_count" ;
		}
	}
}
else {
	foreach($fields AS $item) {
		$field_value_count = $field_value_count + 1;
		if (in_array($item, $field_values)) {
			if ($item == "genecode") { $field_array[] = "geneCode" ;}
			elseif ($item == "laborator") { $field_array[] = "labPerson" ;}
			//elseif ($item == "creation_date") { $field_array[] = "dateCreation" ;}
			else { $field_array[] = $item ;}
		}
		else {
				$errorList[] = "Invalid field header <b>$item</b> in column <b>$field_value_count</b>" ;
		}
	}
}
unset($item);

// check for errors
// if none found ...
if (sizeof($errorList) != 0 ) {
	$title = "$config_sitename: Add Records Error";
	// print html headers
	include_once('../includes/header.php');
	//print errors
	show_errors($errorList);

}
else {
	//checking that each row have the right amount of columns
	foreach($lines AS $item) {
		$x = $x+1;
		$item_columns = explode("	", $item);
		$num_columns = count($item_columns);
		if ($field_value_count != $num_columns) {
			$errorList[] = "Wrong number of columns in row <b>$x</b>." ;
		}
	}
	unset($item,$item_columns);

	//find what column code is in and checks values for existing codes
	$where_code = array_search("code", $field_array);
	//creating a code array for finding duplicates ($code_array)
	foreach($lines AS $item) {
		$item_columns = explode("	", $item);
		$code_array[] = trim($item_columns[$where_code]);
	}
	unset($item,$item_columns);

	if ($input_type == 'seq') {
		//find what column geneCode is in and checks values for existing geneCodes
		$where_geneCode = array_search("geneCode", $field_array);
		//creating a geneCode array for finding if proposed genes exist in db
		foreach($lines AS $item) {
			$item_columns = explode("	", $item);
			$geneCode_array[] = trim($item_columns[$where_geneCode]);
		}
		unset($item,$item_columns);
	}

	// checks for duplicate codes in the list for voucher data
	if ($input_type == 'vouch') {
		$code_list = $code_array;
		foreach($code_list AS $item) {
			$act_code = array_shift($code_list);
			$found_dupl_codes = array_keys($code_list, $item);
			$dupl_codes_act_pos = array_keys($code_array, $item);
			if ($found_dupl_codes) {
				$code_dup_error = "There is duplicate codes in the dataset - can't proceed without unique codes!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Code: <b>$item</b> is found in rows: " ;
				foreach ($dupl_codes_act_pos as $finding) {
					$finding_row = $finding +1;
					$code_dup_error .= "$finding_row ";
				}
				$errorList[] = $code_dup_error;
			}
			foreach ($found_dupl_codes as $code_to_del){
				unset($code_list[$code_to_del]);
			}
		}
		unset($item);
	}
	else { // checks that geneCodes are valid - checks for presence and compares to geneCode's in genes table
		$xrow_num = 0;
		foreach($geneCode_array AS $item) {
			$xrow_num = $xrow_num + 1;
			if ($item == '') {
				$errorList[] = "Invalid entry: <b>Gene code = $item</b> in $xrow_num";
			}
			else {
				@$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
				//select database
				mysql_select_db($db) or die ('Unable to content');
				if( function_exists(mysql_set_charset) ) {
					mysql_set_charset("utf8");
				}
				$querygC = "SELECT geneCode FROM ". $p_ . "genes WHERE geneCode='$item'";
				$resultgC = mysql_query($querygC) or die ("Error in query: $querygC. " . mysql_error());
				//check for empty edits field
				if (mysql_num_rows($resultgC) == 0) {
					$errorList[] = "Gene Code: <b>$item</b> in row: $xrow_num do not exist in gene table
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please add the gene first";
				}
			}
		}
		unset($item, $xrow_num);
		//Checking for duplicate code+Genecode pairs and gives error message if found
		$codeandgeneCode_array = $lines;
		$codeandgeneCode_test = array();
		foreach($codeandgeneCode_array AS $item) {
			$item_columns = explode("	", $item);
			$item_cols = array();
			$item_cols[] = strtolower($item_columns[$where_code]);
			$item_cols[] = strtolower($item_columns[$where_geneCode]);
			$codeandgeneCode_test[] = implode(" + ", $item_cols);
		}
		unset ($dup_c_gc, $item, $item_columns, $item_cols, $xrow_num, $yrow_num);
		$codeandgeneCode_List = $codeandgeneCode_test;
		foreach($codeandgeneCode_test AS $item) {
			$act_code = array_shift($codeandgeneCode_test);
			$found_dupl_codes = array_keys($codeandgeneCode_test, $item);
			$dupl_codes_act_pos = array_keys($codeandgeneCode_List, $item);
			if ($found_dupl_codes) {
				$code_dup_error = "There is duplicate code+genecode pairs in the dataset - can't proceed without unique such pairs!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Code + Genecode: <b>$item</b> is found in rows: " ;
				foreach ($dupl_codes_act_pos as $finding) {
					$finding_row = $finding +1;
					$code_dup_error .= "$finding_row ";
				}
				$errorList[] = $code_dup_error;
			}
		}
		unset($item);
	}

	// checks for existing codes for updating - for seq code+genecode combo
	$xrow_num = 0;
	foreach($lines AS $item) {
		$xrow_num = $xrow_num + 1;
		$item_columns = explode("	", $item);
		if (trim($item_columns[$where_code]) == ''){
			$errorList[] = "Invalid entry: <b>Code</b> in row <b>$xrow_num.</b>";
		}
		else {
			// open database connection
			@$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
			//select database
			mysql_select_db($db) or die ('Unable to content');
			if( function_exists(mysql_set_charset) ) {
				mysql_set_charset("utf8");
			}
			if ($input_type == "vouch") {
				$query_code_test = "SELECT code FROM ". $p_ . "vouchers WHERE code='$item_columns[$where_code]'";
			}
			else {
				$query_code_test = "SELECT code FROM ". $p_ . "sequences WHERE code='$item_columns[$where_code]' AND geneCode='$item_columns[$where_geneCode]'";
			}
			$result_code_test = mysql_query($query_code_test) or die("Error in query: $query. " . mysql_error());
			// if records present
			$num_rows_code_test = mysql_num_rows($result_code_test);
			if( $num_rows_code_test == 0 ) {
				if ($input_type == "vouch") {
					$errorList[] = "Code $item_columns[$where_code] in row <b>$xrow_num</b> doesn't exist!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Can't add values to a non-existing voucher!</b>" ;
				}
				else {
					$errorList[] = "A sequence for code <b>$item_columns[$where_code]</b> and <b>$item_columns[$where_geneCode]</b> 
									in row <b>$xrow_num</b> doesn't exist!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Can't add values to a non-existing sequence!</b>" ;
				}
			}
		}
	}
	unset($item, $xrow_num,$item_columns);

//validate text input fields code, longitude and latitude
if ($input_type == 'vouch') {
	// checking extraction date field field for non-date formatted values and in that case generating error
	$where_extrDate = array_search("dateExtraction", $field_array);
	$xrow_num = 0;
	if ($where_extrDate != False){
		foreach($lines AS $item) {
			$code_num = $xrow_num;
			$xrow_num = $xrow_num + 1;
			$item_columns = explode("	", $item);
			$extrDate = trim($item_columns[$where_extrDate]);
			$dElist = explode("-", $extrDate);
			if ($extrDate != '' && checkdate($dElist[1],$dElist[2],$dElist[0]) == FALSE){
				$errorList[] = "Invalid entry: <b>Extraction Date</b> in row <b>$xrow_num</b>(code = $code_array[$code_num])
					</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Date format should be: YYYY-MM-DD";
			}
		}
	}
	unset($item, $xrow_num,$item_columns);
	// checking latitude field for ilegal values and in that case generating error
	//$where_lat = array_search("latitude", $field_array);
	$where_lat = array_search("latitude", $field_array);
	if ($where_lat) {
		$xrow_num = 0;
		foreach($lines AS $item) {
			$code_num = $xrow_num;
			$xrow_num = $xrow_num + 1;
			$item_columns = explode("	", $item);

			$latitude = $item_columns[$where_lat];
			$latitude = strtoupper($latitude); // to upper in case it is "null"

			// Latitude can be empty space or NULL when user does not have the data. If users are forced to introduce any number
			// we will end up showing maps with spurious coordinates
			// Allow for empty or NULL coordinates
			if( trim($latitude) == "" ) {
				// latitude is empty, allow this
				$latitude = "empty";
			}
			if( trim($latitude) == "NULL" ) {
				// latitude is NULL, allow this
				$latitude = "NULL";
			}

			if( $latitude != "empty" && $latitude != "NULL" && $latitude != NULL && !validate_lat($latitude) ) {
				$errorList[] = "Invalid entry: <b>Latitude</b> in row  <b>$xrow_num </b>(code = $code_array[$code_num]) = $item_columns[$where_lat]
				</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Enter latitude coordinates in <i><b>decimal numbers</b></i> with up to 6 decimal digits.
				</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example:<br />13&deg;08'N = <b>13.133333</b>
				</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;13&deg;08'S = <b>-13.133333</b>";
			}
		}
	}
	unset($item, $xrow_num,$item_columns);

	// checking longitude field for ilegal values and in that case generating error
	$where_long = array_search("longitude", $field_array);
	if ($where_long) {
		foreach($lines AS $item) {
			$code_num = $xrow_num;
			$xrow_num = $xrow_num + 1;
			$item_columns = explode("	", $item);

			$longitude = $item_columns[$where_long];
			$longitude = strtoupper($longitude); // to upper in case it is "null"

			// Longitude can be empty space or NULL when user does not have the data. If users are forced to introduce any number
			// we will end up showing maps with spurious coordinates
			// Allow for empty or NULL coordinates
			if( trim($longitude) == "" ) {
				// longitude is empty, allow this
				$longitude = "empty";
			}
			if( trim($longitude) == "NULL" ) {
				// longitude is NULL, allow this
				$longitude = "NULL";
			}

			if( $longitude != "empty" && $longitude != "NULL" && $longitude != NULL && !validate_long($longitude) ) {
				//echo $valid_long;
				$errorList[] = "Invalid entry: <b>Longitude</b> in row <b>$xrow_num </b>(code = $code_array[$code_num]) = $item_columns[$where_long]
				</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Enter longitude coordinates in <i><b>decimal numbers</b></i> with up to 6 decimal digits.
				</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example:<br />13&deg;08'N = <b>13.133333</b>
				</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;13&deg;08'S = <b>-13.133333</b>";
			}
		}
	}
	unset($item, $xrow_num,$item_columns);
}
else {
	// checking sequences field for ilegal values and in that case generating error
	$where_seq = array_search("sequences", $field_array);
	
	foreach ($lines AS $item){
		$xrow_num = $xrow_num + 1;
		$item_columns = explode("	", $item);
		if (trim($item_columns[$where_seq]) == ''){
			$errorList[] = "Invalid entry: <b>Sequence</b> in row: $xrow_num for code: $item_columns[$where_code] and gene: $item_columns[$where_geneCode].";
		}
	}
	unset($item, $xrow_num,$item_columns);
}
//for updating, checking all values against db to see if some field is already filled

foreach ($field_array as $upfield) {
	if ($upfield !== "code" && $upfield !== "geneCode") {
		$where_up = array_search($upfield, $field_array);
		$xrow_num = 0;
		foreach($lines AS $item) {
			$xrow_num = $xrow_num + 1;
			$item_columns = explode("	", $item);
			if (trim($item_columns[$where_up]) !== '') {
				// open database connection
				@$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
				//select database
				mysql_select_db($db) or die ('Unable to content');
				if( function_exists(mysql_set_charset) ) {
					mysql_set_charset("utf8");
				}
				if ($input_type == 'vouch') { $which_table = "vouchers";}
				else {
					if ($upfield == "primer1" || $upfield == "primer2" || $upfield == "primer3" || $upfield == "primer4"|| $upfield == "primer5" || $upfield == "primer6" ) {
						$which_table = "primers";
					}
					else { $which_table = "sequences";}
				}
				if ($input_type == 'seq') {
					$query_code_test = "SELECT $upfield FROM ". $p_ . "$which_table WHERE code='$item_columns[$where_code]' AND geneCode='$item_columns[$where_geneCode]'";
				}
				else {
					$query_code_test = "SELECT $upfield FROM ". $p_ . "$which_table WHERE code='$item_columns[$where_code]'";
				}
				$result_code_test = mysql_query($query_code_test) or die("Error in query: $query_code_test. " . mysql_error());
				// if records present
				$num_rows_code_test = mysql_num_rows($result_code_test);
				if( $num_rows_code_test > 0 ) {
					while ($row = mysql_fetch_object($result_code_test)) {
						if (trim($row->$upfield) !== ''){
							if ($input_type == "vouch"){
								$errorList[] = "Code <b>$item_columns[$where_code]</b> in row <b>$xrow_num</b> already have information in field <b>$upfield</b>!
												</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
												Please edit this manually!</b>" ;
							}
							else {
								$errorList[] = "Code <b>$item_columns[$where_code]</b> for gene <b>$item_columns[$where_geneCode]</b> in row <b>$xrow_num</b> already have information in field <b>$upfield</b>!
												</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
												Please edit this manually!</b>" ;
							}
						}
					}
				}
			}
		}
	}
}


//check for error and if none proceed with inputting into mysql db
if (sizeof($errorList) != 0 ) {
	$title = "$config_sitename: Add Records Error";
	// print html headers
	include_once('../includes/header.php');
	//print errors
	show_errors($errorList);
}
else {
	foreach($lines AS $item) {
		unset($prexists, $qadd_primers, $qup_primers, $qup_fields, $qup_primer_fields, $q_values, $q_fields, $q_primers, $fields_to_query, $created_record, $q_primer_fields, $q_primer_values, $q_field_array, $if_primers, $q_edit_add);
		// open database connection
		@$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
		//select database
		mysql_select_db($db) or die ('Unable to content');
		if( function_exists(mysql_set_charset) ) {
			mysql_set_charset("utf8");
		}


		if( $php_version == "5" ) {
			date_default_timezone_set($date_timezone); //php5
		}

		//setting the edits add values
		$latesteditor = utf8_encode($_SESSION['SESS_FIRST_NAME']. " ". $_SESSION['SESS_LAST_NAME']);
		$edit_columns = explode("	", $item);
		if ($input_type == 'vouch') {
			$editsadd = "Voucher information updated by ". $latesteditor ." on ";
		}
		else {
			if (array_search("primer1", $field_array) || array_search("primer2", $field_array) ||array_search("primer3", $field_array) ||array_search("primer4", $field_array) ||array_search("primer5", $field_array) || array_search("primer6", $field_array) ) {
				$editsadd = $edit_columns[$where_geneCode] . " sequence and primer information updated by " . $latesteditor . " on ";
			}
			else {
				$editsadd = $edit_columns[$where_geneCode] . " sequence information updated by " . $latesteditor . " on ";
			}
		}

		//mysql_query("time for add-list");
		$querytime = "SELECT NOW()";
		$resulttime = mysql_query($querytime) or die ("Error in query: $querytime. " . mysql_error());
		$rowtime    = mysql_result($resulttime,0);
		$editsadd = $editsadd . $rowtime . "\n";
		$queryed = "SELECT edits FROM ". $p_ . "vouchers WHERE code='$edit_columns[$where_code]'";
		$resulted = mysql_query($queryed) or die ("Error in query: $queryed. " . mysql_error());
		//check for empty edits field
		if (mysql_num_rows($resulted) > 0) {
			$rowed    = mysql_result($resulted,0);
			//adding to edits
			$editsadd = $editsadd . $rowed ;
		}
			
		// generate and execute query
		if ($input_type == 'vouch') {
			$qup_start = "UPDATE ". $p_ . "vouchers SET ";
		}
		else {
			$qup_start = "UPDATE ". $p_ . "sequences SET "; 
			if (array_search("primer1", $field_array) || array_search("primer2", $field_array) ||array_search("primer3", $field_array) ||array_search("primer4", $field_array)) {
				$querypr = "SELECT id FROM ". $p_ . "primers WHERE code='$edit_columns[$where_code]' AND geneCode='$edit_columns[$where_geneCode]'";
				$resultpr = mysql_query($querypr) or die ("Error in query: $querypr. " . mysql_error());
				if (mysql_num_rows($resultpr) > 0) {
					$prexists = "TRUE";
				}
				else {
					$prexists = "FALSE";
				}
			}
		}
		//$qup_fields = array();$qup_primer_fields = array();
		foreach ($field_array as $field) {
			$item_columns = explode("	", $item); // enter field value into query
			$where_field = array_search($field, $field_array);
			if (trim($item_columns[$where_field]) !== '') {
				if ($field == "code") { 
					$created_record[] = $item_columns[$where_field];
					$actual_code = $item_columns[$where_field];
				} // just adds the code to a list to show as successful upon finished creations, skips adding code to update query
				elseif ($field == "geneCode") { 
					$created_record[] = $item_columns[$where_field];
					$actual_geneCode = $item_columns[$where_field];
				} // just adds the genecode to a list to show as successful upon finished creations, skips adding genecode to update query
				else {
					if ($field == 'sequences') {
						$cleaned_sequences = preg_replace("/\s/", "",trim($item_columns[$where_field]));
						$cleaned_sequences = strtoupper($cleaned_sequences);
						$qup_fields[] = "$field='$cleaned_sequences'";
					}
					elseif ($field == 'primer1' || $field == 'primer2' || $field == 'primer3' || $field == 'primer4' || $field == 'primer5' || $field == 'primer6') {
						$if_primers = $if_primers + 1;
						$q_primer_fields[] = $field;
						$q_primer_values[] = "'$item_columns[$where_field]'";
						$qup_primer_fields[] = "$field='$item_columns[$where_field]'";
					}
					else {
						// check is value is empty
						if( trim($item_columns[$where_field]) == "" ) {
							$q_value = "NULL";
						}
						// check is value is "null"
						elseif( trim(strtoupper($item_columns[$where_field])) == "NULL" ) {
							$q_value = "NULL";
						}
						else {
							$q_value = "$item_columns[$where_field]";
						}
						$qup_fields[] = "$field='$q_value'";
					}
				}
			}
		}
		unset($field);

		$created_records[] = implode(" ", $created_record);
		if ($if_primers > 0) {
			if (isset($q_primer_fields)) {$q_primer_fields = implode(", " , $q_primer_fields);}
			if (isset($q_primer_values)) {$q_primer_values = implode(", " , $q_primer_values);}
			if (isset($qup_primer_fields)) {$qup_primer_fields = implode(", " , $qup_primer_fields);}
		}
		if (isset($qup_fields)) {$qup_fields = implode (", " , $qup_fields);}

		if ($input_type == 'vouch') {
			if (isset($qup_fields)) {
				$q_fields = "$qup_start $qup_fields, timestamp=NOW(), edits='$editsadd', latesteditor='$latesteditor' WHERE code='$actual_code'";
			}
		}
		else {
			$dateModification = date('Y-m-d');
			if (isset($qup_fields)) {
				$q_fields = "$qup_start $qup_fields, timestamp=NOW(), dateModification='$dateModification' WHERE code='$actual_code' AND geneCode='$actual_geneCode'";
			}
			if ($if_primers > 0) {
				if ($prexists == "FALSE") {
					if (isset($q_primer_fields) && isset($q_primer_values)){
						$qadd_primers .= "INSERT INTO ". $p_ . "primers ($q_primer_fields, code, geneCode, timestamp) VALUES ( $q_primer_values, '$item_columns[$where_code]', '$item_columns[$where_geneCode]', NOW())";
						$q_primers = $qadd_primers;
					}
				}
				else {
					if (isset($qup_primer_fields)){
						$qup_primers .= "UPDATE ". $p_ . "primers SET $qup_primer_fields, timestamp=NOW() WHERE code='$actual_code' AND geneCode='$actual_geneCode'";
						$q_primers = $qup_primers;
					}
				}
			}
		}
		
		//echo "$q_fields<br />"; //only for test purposes
		//echo "$q_primers<br />"; //only for test purposes
		if (isset($q_fields)){
			$result = mysql_query($q_fields) or die ("Error in query: $q_fields. " . mysql_error());
		}
		if ($input_type == 'seq' && isset($q_fields) || isset($q_primers)) {
				$q_edit_add = "UPDATE ". $p_ . "vouchers SET edits='$editsadd', latesteditor='$latesteditor', timestamp=NOW() WHERE code='$item_columns[$where_code]'";
				$result = mysql_query($q_edit_add) or die ("Error in query: $q_edit_add. " . mysql_error());
		}
		if ($if_primers > 0 && isset($q_primers)) {
			$result = mysql_query($q_primers) or die ("Error in query: $q_primers. " . mysql_error());
		}
	}

	// process title
	$title = "$config_sitename - Records updated";

	// print html headers
	include_once('../includes/header.php');

	// print navegation bar
	admin_nav();

	// begin HTML page content
	echo "<div id=\"content_narrow\">";
	echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
			<tr><td valign=\"top\">";
	// print result
	$created_records = implode("," , $created_records);
	echo "<span class=\"title\"><img src=\"images/success.png\" alt=\"\"> Record update was successful for the following codes: $created_records.</span>";
	echo "<td class=\"sidebar\" valign=\"top\">";
	admin_make_sidebar();
	echo "</td>";
	echo "</tr>
		</table> <!-- end super table -->
		</div> <!-- end content -->";
	make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);

	echo "</body>";
	echo "</html>";
	

	}
}

?>
