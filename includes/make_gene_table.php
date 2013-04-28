<?php
// #################################################################################
// #################################################################################
// Voseq includes/process_dataset.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Processes info from ../create_dataset.php to create dataset
// #################################################################################
// #################################################################################
// Section: Startup/includes
// #################################################################################
//check admin login session
include'../login/auth.php';
// includes
ob_start();//Hook output buffer - disallows web printing of file info...
include_once'../conf.php';
ob_end_clean();//Clear output buffer//includes
include '../functions.php';
include '../markup-functions.php';
include "translation_functions.php";

// set increased time limit for larger calculations -> 1000 seconds
set_time_limit(1500);

$warning = array();
$charset_count = array();
$errorList = array();
$geneCodes = array();
$positions = array();
$rfs = array();

// open database connections
@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect');
mysql_select_db($db) or die ('Unable to select database');
if( function_exists(mysql_set_charset)) {
	mysql_set_charset("utf8");
}
// #################################################################################
// Section: Functions - clean_item() and show_errors()
// #################################################################################
function clean_item ($item) {
	$item = stripslashes($item);
	$item = str_replace("'", "", $item);
	$item = str_replace('"', "", $item);
	$item = str_replace(',', "", $item);
	$item = preg_replace('/^\s+/', '', $item);
	$item = preg_replace('/^\t+/', '', $item);
	$item = preg_replace('/\s+$/', '', $item);
	$item = strtolower($item);
	return $item;
}

function show_errors($se_in) {
		// error found
			// print navegation bar
			nav();
			// begin HTML page content
			echo "<div id=\"content_narrow\">";
			echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
					<tr><td valign=\"top\">";
		// print as list
		echo "<img src=\"../images/warning.png\" alt=\"\"> The following errors were encountered:";
		echo '<br>';
		echo '<ul>';
			$se_in = array_unique($se_in);
			$se_in[] = "</br>Please revise your data!"; 
		foreach($se_in AS $item) {
			echo "$item</br>";
		}
		echo "</td>";
		
		echo "<td class='sidebar'>";
		make_sidebar(); 
		echo "</td>";
		echo "</tr>
				</table> <!-- end super table -->
				</div> <!-- end content -->";
		//make footer
		make_footer($date_timezone, $config_sitename, $version, $base_url);
		?></body></html><?php
}
function nucfreqs($sequence) { // gives array with numbers for occurrences of each nucleotide as well as ? - N etc
	$sequence_array = preg_split('#(?<=.)(?=.)#s', $sequence);
	$freqvalues = array_count_values($sequence_array);
	return $freqvalues;
}
function nuc_chars($nuc_count_array){ // calculates numbers for a nuc_count array
	$nchars = $nuc_count_array['A'] + $nuc_count_array['C'] + $nuc_count_array['T'] + $nuc_count_array['G'];
	$totalchars = array_sum($nuc_count_array);
	$dcomp = $nchars/$totalchars;
	//$table1[$geneCode][] = $dcomp;
	$missingchars = $totalchars - $nchars;
	return array($totalchars, $nchars, $missingchars);
}
function PIchars ($sequence_arr) { // gives frequencies of variable, PI and conserved chars for a set of sequences
	$length_seq = count(preg_split('#(?<=.)(?=.)#s', $sequence_arr[0]));
	$V = $C = $PI = 0;
	for ($i=0;$i <= $length_seq; $i++){
		$prfreqs = $cleaned_pfreqs = $cleaned_pfreqs2 = $posarr = array();
		foreach ($sequence_arr as $ind_seq) {
			$sarr = preg_split('#(?<=.)(?=.)#s', $ind_seq); // making sequence/nucleotide array
			$posarr[] = $ind_seq[$i];
		}
		$pfreqs = array_count_values($posarr);
		foreach (array('A', 'T', 'C', 'G') as $base){
			if ($pfreqs[$base] > 0){ 
				$cleaned_pfreqs[$base] = $pfreqs[$base];
			}
			if ($pfreqs[$base] > 1){ 
				$cleaned_pfreqs2[$base] = $pfreqs[$base];
			}
		}
		$cleaned_length = array_sum($cleaned_pfreqs);
		if (count($cleaned_pfreqs) < 2) {$C = $C+1;}
		else{
			if (count($cleaned_pfreqs2) > 1) { $PI=$PI+1;$V=$V+1;}
			else{ $V=$V+1;}
		
		}
	}
	$V = $V/$length_seq;
	$C = $C/$length_seq;
	$PI = $PI/$length_seq;
	return array($V, $PI, $C );
}
// #################################################################################
// Section: Get code(s) and gene(s)
// #################################################################################

// #################################################################################
// Section: Set variables
// #################################################################################
if (isset($_POST['geneCodes'])){
	foreach ( $_POST['geneCodes'] as $k1=> $c1){ //putting choosen genes into array
		if ($c1 == 'on')	{
			$geneCodes[] =  $k1;
		}
	}
}
elseif ($_POST['genesets'] == "Choose geneset"){
	$errorList[] = "No genes choosen - Please try again!"; 
}
// checking geneset choice
$geneset = $_POST['genesets'];
$geneset_taxa = array();
if ($geneset != "Choose geneset"){
	$TSquery = "SELECT geneset_list FROM ". $p_ . "genesets WHERE geneset_name='$geneset'";
	$TSresult = mysql_query($TSquery) or die("Error in query: $TSquery. " . mysql_error());
		// if records present
		
		if( mysql_num_rows($TSresult) > 0 ) {
			while( $TSrow = mysql_fetch_object($TSresult) ) {
				$geneset_taxa = explode(",", $TSrow->geneset_list );
			}
		}
	else {$errorList[] = "No gene set named <b>$geneset</b> exists in database!";}
}else {unset($geneset_taxa);}

// merging choosen gene set taxa and input taxa lists
if (isset($geneset_taxa) && isset($geneCodes)){$geneCodes = array_merge( $geneset_taxa, $geneCodes) ;}
elseif (isset($geneset_taxa) && ! isset($geneCodes)){$geneCodes = $geneset_taxa ;}
elseif (! isset($geneset_taxa) && isset($geneCodes)){$geneCodes = $geneCodes ;}
else { $errorList[] = "No genes are chosen!</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Pointless to make a table without genes..."; }


 //input data
if (trim($_POST['codes']) != "") {
	$raw_codes = explode("\n", $_POST['codes']);
}
else { 
	unset($raw_codes); 
}
if ( $_POST['field_delimitor'] == 'comma') {
	$field_delimitor = ",";
	}
else {
	$field_delimitor = "	";
}
if ( $_POST['decimal'] == 'comma') {
	$decimal = ",";
	}
else {
	$decimal = ".";
}
if ( $_POST['codpos'] == 'yes') {
	$codpos = "yes";
	}
else {
	$codpos = "no";
}
$result = mysql_query("set names utf8") or die("Error in query: $query. " . mysql_error());

$seq_string = "?"; #this will be a replacement for NULL sequences
$bp = 0; #number of base pairs


foreach ($geneCodes AS $item) { // building full dataset bp count + setting charset_count[] values and setting reading frames
	//$item = clean_item($item);
	$gCquery = "SELECT geneCode, length, readingframe FROM ". $p_ . "genes WHERE geneCode='$item'";
	$gCresult = mysql_query($gCquery) or die("Error in query: $query. " . mysql_error());
	// if records present
	if( mysql_num_rows($gCresult) > 0 ) {
		while( $row = mysql_fetch_object($gCresult) ) {
			if ($row->length != "0") {
				$bp = $bp + $row->length;
				$charset_count[$item] = $row->length;
			}
			else { 
				if ($format != "FASTA"){
					$errorList[] = "The length of gene <b>$item</b> (in bp's) has not been specified!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Please do that in the gene edit section!";
				}
			}
			$rf = $row->readingframe;
			$rfs[$item] = $rf ;
		}
	}
	else {
		$errorList[] = "Gene <b>$item</b> does not exist in database!</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please add it in the gene edit section!";
	}
}
unset($item);

// checking taxonset choice
$taxonset = $_POST['taxonsets'];
$taxonset_taxa = array();
if ($taxonset != "Choose taxonset"){
	$TSquery = "SELECT taxonset_list FROM ". $p_ . "taxonsets WHERE taxonset_name='$taxonset'";
	$TSresult = mysql_query($TSquery) or die("Error in query: $TSquery. " . mysql_error());
		// if records present
		
		if( mysql_num_rows($TSresult) > 0 ) {
			while( $TSrow = mysql_fetch_object($TSresult) ) {
				$taxonset_taxa = explode(",", $TSrow->taxonset_list );
			}
		}
	else {$errorList[] = "No taxon set named <b>$taxonset</b> exists in database!";}
}else {unset($taxonset_taxa);}

// merging choosen taxon set taxa and input taxa lists
if (isset($taxonset_taxa) && isset($raw_codes)){$raw_codes = array_merge( $taxonset_taxa, $raw_codes) ;}
elseif (isset($taxonset_taxa) && ! isset($raw_codes)){$raw_codes = $taxonset_taxa ;}
elseif (! isset($taxonset_taxa) && isset($raw_codes)){$raw_codes = $raw_codes ;}
else { $errorList[] = "No taxa are chosen!</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Pointless to make a table without taxa..."; }



$codes = array();
if (isset($raw_codes)){
$raw_codes = array_unique($raw_codes); 
foreach($raw_codes AS $item) {
	if ( clean_item($item) != "") {
		$item = clean_item($item);
		$item = trim($item);
		if (strpos($item, "--") === 0) {$item = str_replace("--","",$item);}
		$cquery = "SELECT code FROM ". $p_ . "vouchers WHERE code='$item'";
		$cresult = mysql_query($cquery) or die("Error in query: $query. " . mysql_error());
		// if records present
		if( mysql_num_rows($cresult) > 0 ) {
			while( $row = mysql_fetch_object($cresult) ) {		
				array_push($codes, $item);
			}
		}
		else {
		$errorList[] = "No voucher named <b>$item</b> exists in database!</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please add it in the voucher section or remove it from taxon set!";
		}
	}
}unset($item);
$codes = array_unique($codes);
}

// removing taxa from list that have the removal code -- before them
$codes_to_remove = array();
if (isset($taxonset_taxa) && isset($raw_codes)){
	$raw_codes_delete = $raw_codes;
	foreach($raw_codes_delete AS $item) {
		if(strpos($item,'--') !== false) {
			$item = clean_item($item);
			$item = trim($item);
			$item2 = str_replace('--','',$item);
			$codes_to_remove[] = $item2;
		}
	}
	$codes = array_diff($codes, $codes_to_remove);
}unset($item,$item2);



//if no taxa left or presented
$number_of_taxa = count($codes);
if ($number_of_taxa == 0) {$errorList[] = "No codes specified! No use creating empty datasets...
											</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
											Please go back and add voucher codes to run!"; 
}
// check for gene and taxon sets missing all sequences
foreach ($geneCodes AS $geneCode){
	$test = 0;
	foreach ($codes as $c){
		$cquery = "SELECT sequences FROM ". $p_ . "sequences WHERE code='$c' AND geneCode='$geneCode'";
		$cresult = mysql_query($cquery) or die("Error in query: $query. " . mysql_error());
		// if records present
		if( mysql_num_rows($cresult) > 0 ) {
			while( $row = mysql_fetch_object($cresult) ) {
				$cseq = $row->sequences;
					$cseq = str_replace(array("?", "-"), "", $cseq);
					if (strlen($cseq) > 0 ) {$test = 1;}
			}
		}
	}
	if ($test == 0){ 
		$errorList[] = "Gene $geneCode has no sequences for choosen taxa!
						</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
						Cannot calculate values for it!"; 
	}
}
// #################################################################################
// Section: check for error and if none proceed with building dataset
// #################################################################################
if (sizeof($errorList) != 0 ){
	$title = "$config_sitename: Dataset Error";
	// print html headers
	$admin = false;
	$in_includes = true;
	include_once 'header.php';
	//print errors
	show_errors($errorList);
}
else{ 
// #################################################################################
// Section: Start building dataset for calculations
// #################################################################################

// #################################################################################
	// Section: Build dataset and adding choosen name-extensions
	// #################################################################################

	$num_genes = 0;
	$output_lines = "";
	$geneinfo_array = array();
	$taxout_array = array();
	$seqout_array = array();
	foreach ($geneCodes AS $geneCode) {
		$num_genes = $num_genes + 1;
		$query_b = "SELECT genetype, prot_code, aligned, length, intron, readingframe FROM ". $p_ . "genes WHERE geneCode='$geneCode'";
		$result_b = mysql_query($query_b) or die("Error in query: $query_b. " . mysql_error());
		// if records present
		if( mysql_num_rows($result_b) > 0 ) {
			while( $generow = mysql_fetch_object($result_b) ) {
			$geneinfo_array[$geneCode]['genetype'] = $generow->genetype;
			$geneinfo_array[$geneCode]['aligned'] = $generow->aligned;
			$geneinfo_array[$geneCode]['length'] = $generow->length;
			$geneinfo_array[$geneCode]['introns'] = $generow->intron;
			$geneinfo_array[$geneCode]['reading_frame'] = $generow->readingframe;
			$geneinfo_array[$geneCode]['prot_code'] = $generow->prot_code;
			}
		}
		else {$errorList[] = "couldnt fetch result..."; }
		
		// If intron - get new length without it
		if (isset($geneinfo_array[$geneCode]['introns']) && $geneinfo_array[$geneCode]['introns'] != '' && $geneinfo_array[$geneCode]['introns'] != 'NULL'){
			$geneinfo_array[$geneCode]['intron_info'] = remove_introns(str_pad("A", $charset_count[$geneCode], "A"),$geneinfo_array[$geneCode]['introns'] );
			$geneinfo_array[$geneCode]['length'] = $geneinfo_array[$geneCode]['intron_info'][1];
		}
		if ($geneinfo_array[$geneCode]['prot_code'] == "yes" && $geneinfo_array[$geneCode]['reading_frame'] == "" && $codpos == 'yes'){
				$errorList[] = "The $geneCode is set as Protein Coding but doesnt have a specified reading frame!
							</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please edit gene info.";
		}
		if ($geneinfo_array[$geneCode]['aligned'] == 'yes' ){ // checking if aligned or not
			foreach ($codes AS $item) {
				$item = clean_item($item);
				$taxout_array[$geneCode][$item] = $item;
				// #################################################################################
				// Section: Sequence builder
				// #################################################################################
									
				$query_b = "SELECT sequences FROM ". $p_ . "sequences WHERE code='$item' AND geneCode='$geneCode'";
				$result_b = mysql_query($query_b) or die("Error in query: $query_b. " . mysql_error());
				// if records present
				if( mysql_num_rows($result_b) > 0 ) {
					while( $row_b = mysql_fetch_object($result_b) ) {
						$seq = $row_b->sequences;
						if (strlen($charset_count[$geneCode] < strlen($seq))){ // checking for too long sequence
							$errorList[] = "The $geneCode sequence of $item is longer (". strlen($seq) . ">" . $charset_count[$geneCode] .")that the specified gene length!
							</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please edit gene length or check the sequence";
						}
						$seqout_array[$geneCode][$item] = str_pad($seq, $charset_count[$geneCode], "?");
						// if introns capture "cleaned" sequence and get intron info in new table [genecode]['intron_info']
						if (isset($geneinfo_array[$geneCode]['introns']) && $geneinfo_array[$geneCode]['introns'] != ''&& $geneinfo_array[$geneCode]['introns'] != 'NULL'){
							$geneinfo_array[$geneCode]['intron_info'] = remove_introns($seqout_array[$geneCode][$item], $geneinfo_array[$geneCode]['introns']);
							$seqout_array[$geneCode][$item] = $geneinfo_array[$geneCode]['intron_info'][0];
						}
					}
				}
				else {
						//echo "no sequence $geneCode + $item  </br>";
						$seq = "?"; 
						if (isset($geneinfo_array[$geneCode]['intron_info'])){
							$seqout_array[$geneCode][$item] = str_pad($seq, $geneinfo_array[$geneCode]['length'], "?");
						}
						else {
							$seqout_array[$geneCode][$item] = str_pad($seq, $charset_count[$geneCode], "?");
						}
				}
			}
		}
		else {
			if (count($codes) > 1){
				$errorList[] = "You have chosen an unaligned alignment/gene ($geneCode) and more than one taxa!
							</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Cannot make calculations on that!
							</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please edit gene alignment status or number of taxa";
			}
		}
	}
	unset($item);

	// #################################################################################
	// Section: check for error and if none proceed with building dataset
	// #################################################################################
	if (sizeof($errorList) != 0 ){
		$title = "$config_sitename: Dataset Error";
		// print html headers
		$admin = false;
		$in_includes = true;
		include_once 'header.php';
		//print errors
		show_errors($errorList);
	}
	else{ //start calculating and building table
		foreach ($geneCodes AS $geneCode) {
			//  build position lists
			foreach ($codes as $item) {
				$geneinfo_array[$geneCode]['tot'][] =  $seqout_array[$geneCode][$item];
				
				if ($geneinfo_array[$geneCode]['prot_code'] == "yes" && $codpos == 'yes'){
					$sequence_array = preg_split('#(?<=.)(?=.)#s', $seqout_array[$geneCode][$item]); // making sequence/nucleotide array
					for ($i=0; $i <= $geneinfo_array[$geneCode]['length']; $i+=3){
						if ($geneinfo_array[$geneCode]['reading_frame'] == "1") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['1'][] = $sequence_array[$i];
							}
						}
						if ($geneinfo_array[$geneCode]['reading_frame'] == "2") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['3'][] = $sequence_array[$i];
							}
						}
						if ($geneinfo_array[$geneCode]['reading_frame'] == "3") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['2'][] = $sequence_array[$i];
							}
						}
					}
					for ($i=1; $i <= $geneinfo_array[$geneCode]['length']; $i+=3){
						if ($geneinfo_array[$geneCode]['reading_frame'] == "1") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['2'][] = $sequence_array[$i];
							}
						}
						if ($geneinfo_array[$geneCode]['reading_frame'] == "2") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['1'][] = $sequence_array[$i];
							}
						}
						if ($geneinfo_array[$geneCode]['reading_frame'] == "3") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['3'][] = $sequence_array[$i];
							}
						}
					}
					for ($i=2; $i <= $geneinfo_array[$geneCode]['length']; $i+=3){
						if ($geneinfo_array[$geneCode]['reading_frame'] == "1") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['3'][] = $sequence_array[$i];
							}
						}
						if ($geneinfo_array[$geneCode]['reading_frame'] == "2") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['2'][] = $sequence_array[$i];
							}
						}
						if ($geneinfo_array[$geneCode]['reading_frame'] == "3") {
							if (isset($sequence_array[$i])){
								$geneinfo_array[$geneCode][$item]['1'][] = $sequence_array[$i];
							}
						}
					}
					$geneinfo_array[$geneCode]['1'][] = implode("",$geneinfo_array[$geneCode][$item]['1']);
					$geneinfo_array[$geneCode]['2'][] = implode("",$geneinfo_array[$geneCode][$item]['2']);
					$geneinfo_array[$geneCode]['3'][] = implode("",$geneinfo_array[$geneCode][$item]['3']);
				}
			}
			// combine to sets
			$geneinfo_array[$geneCode]['totset'] = implode("", $geneinfo_array[$geneCode]['tot']);
			if ($geneinfo_array[$geneCode]['prot_code'] == "yes" && $codpos == 'yes'){
				$geneinfo_array[$geneCode]['1set'] = implode("", $geneinfo_array[$geneCode]['1']);
				$geneinfo_array[$geneCode]['2set'] = implode("", $geneinfo_array[$geneCode]['2']);
				$geneinfo_array[$geneCode]['3set'] = implode("", $geneinfo_array[$geneCode]['3']);
			}			
			// start calculations of base frequencies
			if ($geneinfo_array[$geneCode]['prot_code'] == "yes" && $codpos == 'yes'){
				foreach (array("tot", "1", "2", "3") as $i){
					$iset = $i . "set"; $i_count = $i . "_count";
					//$sequence_array = preg_split('#(?<=.)(?=.)#s', $geneinfo_array[$geneCode][$iset]);
					$geneinfo_array[$geneCode][$i_count] = nucfreqs($geneinfo_array[$geneCode][$iset]);
				}
			}
			else {
				$geneinfo_array[$geneCode]['tot_count'] = nucfreqs($geneinfo_array[$geneCode]['totset']);
				//$sequence_array = preg_split('#(?<=.)(?=.)#s', $geneinfo_array[$geneCode]['totset']);
				//$geneinfo_array[$geneCode]['tot_count'] = array_count_values($sequence_array);
			}
		}
		// Build table
		// Check if any gene is protcoding
		$cpon = "no";
		foreach ($geneCodes as $gcpc){
			if ($geneinfo_array[$gcpc]['prot_code'] == "yes" && $codpos == 'yes'){
				$cpon = "on";
			}
		}
		// Check if any gene has introns
		$intrOn = "no";
		foreach ($geneCodes as $gcpc){
			if ($geneinfo_array[$gcpc]['intron_info'][2] > 0){
				$intrOn = "on";
			}
		}
		// headers
		$table1 = array();
		$table1['headers'][] = "Gene";
		if ($cpon == "on") {$table1['headers'][] = "pos.";}
		$table1['headers'][] = "Gene type";
		$table1['headers'][] = "Length";
		$table1['headers'][] = "Dataset completion (%)";
		
		if (count($codes) > 1){
			$table1['headers'][] = "Variable (%)";
			if (count($codes) > 3){$table1['headers'][] = "Pars. Inf.(%)";}
			$table1['headers'][] = "Conserved (%)";
		}
		$table1['headers'][] = "Freq. A (%)";
		$table1['headers'][] = "Freq. T (%)";
		$table1['headers'][] = "Freq. C (%)";
		$table1['headers'][] = "Freq. G (%)";
		if($intrOn == 'on'){ 
			$table1['headers'][] = "Introns (n)";
			$table1['headers'][] = "Tot. intron length (bp)";
		}
		// input headers in table
		$xls = implode($field_delimitor, $table1['headers']) . "\n";
		// looped values for genes

		foreach ($geneCodes as $geneCode){
			$table1[$geneCode][] = $geneCode;
			if ($cpon == "on") {$table1[$geneCode][] = "";}
			$table1[$geneCode][] = $geneinfo_array[$geneCode]['genetype'];
			$table1[$geneCode][] = $geneinfo_array[$geneCode]['length'];
			$capt_array = nuc_chars($geneinfo_array[$geneCode]['tot_count']);
			$table1[$geneCode][] = number_format(round(($capt_array[1] / $capt_array[0])*100,1), 2, $decimal, '');
			
			if (count($codes) > 1){
				$capt_array_PI = PIchars($geneinfo_array[$geneCode]['tot']);
				$table1[$geneCode][] = number_format(round($capt_array_PI[0]*100, 2), 2, $decimal, '');
				if (count($codes) > 3){$table1[$geneCode][] = number_format(($capt_array_PI[1]*100), 2, $decimal, '');}
				$table1[$geneCode][] = number_format(100-round($capt_array_PI[0]*100, 2), 2, $decimal, '');
			}
			$table1[$geneCode][] = number_format(($geneinfo_array[$geneCode]['tot_count']['A'] / $capt_array[1])*100, 2, $decimal, '');
			$table1[$geneCode][] = number_format(($geneinfo_array[$geneCode]['tot_count']['T'] / $capt_array[1])*100, 2, $decimal, '');
			$table1[$geneCode][] = number_format(($geneinfo_array[$geneCode]['tot_count']['C'] / $capt_array[1])*100, 2, $decimal, '');
			$table1[$geneCode][] = number_format(($geneinfo_array[$geneCode]['tot_count']['G'] / $capt_array[1])*100, 2, $decimal, '');
			if($intrOn == 'on'){ 
				if (isset($geneinfo_array[$geneCode]['intron_info'])){
					$table1[$geneCode][] = $geneinfo_array[$geneCode]['intron_info'][2];
					$table1[$geneCode][] = $geneinfo_array[$geneCode]['intron_info'][3];
					// for ($j = 3; $j <= count($geneinfo_array[$geneCode]['intron_info']);$j++){
						// $table1[$geneCode][] = $geneinfo_array[$geneCode]['intron_info'][$j];
					// }
				}
				else {
					$table1[$geneCode][] = "";
					$table1[$geneCode][] = "";
				}
			}
			//put gene info in table
			$xls .= implode($field_delimitor, $table1[$geneCode]) . "\n";
			// separate lines for codon positions
			if ($geneinfo_array[$geneCode]['prot_code'] == "yes" && $codpos == 'yes'){
				for ($i=1; $i <= 3; $i++){
					$table1[$i] = array();
					$table1[$i][] = "";
					$table1[$i][] = $i;
					$table1[$i][] = "";
					$table1[$i][] = count($geneinfo_array[$geneCode][$codes[0]][$i]);
					$capt_array = nuc_chars($geneinfo_array[$geneCode][$i."_count"]);
					$table1[$i][] = number_format(($capt_array[1] / $capt_array[0])*100, 2, $decimal, '');
					
					if (count($codes) > 1){
						$capt_array_PI = PIchars($geneinfo_array[$geneCode][$i]);
						$table1[$i][] = number_format(round($capt_array_PI[0]*100,2), 2, $decimal, '');
						if (count($codes) > 3){$table1[$i][] = number_format($capt_array_PI[1]*100, 2, $decimal, '');}
						$table1[$i][] = number_format(100-round($capt_array_PI[0]*100, 2), 2, $decimal, '');
					}
					$table1[$i][] = number_format(($geneinfo_array[$geneCode][$i."_count"]['A'] / $capt_array[1])*100, 2, $decimal, '');
					$table1[$i][] = number_format(($geneinfo_array[$geneCode][$i."_count"]['T'] / $capt_array[1])*100, 2, $decimal, '');
					$table1[$i][] = number_format(($geneinfo_array[$geneCode][$i."_count"]['C'] / $capt_array[1])*100, 2, $decimal, '');
					$table1[$i][] = number_format(($geneinfo_array[$geneCode][$i."_count"]['G'] / $capt_array[1])*100, 2, $decimal, '');
					if($intrOn == 'on'){ 
							$table1[$i][] = "";
							$table1[$i][] = "";
					}
					//$table1[$i][] = implode("-", $geneinfo_array[$geneCode][$codes[0]][$i]);
					//put codon position info in table
					$xls .= implode($field_delimitor, $table1[$i]) . "\n";
				}
			}
		}
	}
	// #################################################################################
	// Section: Create downloadable file with table
	// #################################################################################
	# filename for download
	if( $php_version == "5" ) {
		//date_default_timezone_set($date_timezone);php5
		date_default_timezone_set($date_timezone);
	}
	$excel_file = "db_table_" . date('Ymd') . ".xls";
	header("Content-Disposition: attachment; filename=\"$excel_file\"");
	header("Content-Type: application/vnd.ms-excel");
	echo $xls;
}
?>
