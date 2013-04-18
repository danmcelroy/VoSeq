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
			echo "$item</br></br>";
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



// #################################################################################
// Section: Get code(s) and gene(s)
// #################################################################################

// #################################################################################
// Section: Special mode - set variables
// #################################################################################
// checking to see if special mode is enabled, and in that case just copy the values and fix the by-gene values and proceed to building the dataset
if (isset($_POST['geneCodes2']) && isset($_POST['codes2']) && isset($_POST['gene_by_positions2']) && isset($_POST['gene_positions2'])){
	$format = $_POST['format2']; //echo "format = >$format<</br>";
	$geneCodes = explode(",", $_POST['geneCodes2']);
	$codes = explode(",", $_POST['codes2']); // $codes3 = implode(",",$codes);echo "codes: " . $codes3 ."</br>";
	$taxonadds = explode(",",$_POST['taxonadds2']); //$taxonadds3 = implode(",",$taxonadds);echo "taxonadds: " . $taxonadds3 ."</br>";
	$outgroup = $_POST['outgroup2']; //echo "outgroup: $outgroup</br>";
	$positions = array("special");
	$by_positions = "special";
	$number_of_taxa = $_POST['number_of_taxa'];
	$gene_positions2 = $_POST['gene_positions2'];
	$gene_by_positions = $_POST['gene_by_positions2'];
	$ignore_introns = $_POST['ignore_introns'];
	// $aligned = $_POST['aligned'];
	// $prot_code = $_POST['prot_code'];
	// $genetic_code = $_POST['genetic_code'];
	$special = 'yes'; // for recognization of special mode read
	foreach ($geneCodes as $genecode2){
		//get charset counts and rfs
		$charset_count[$genecode2] = $_POST[$genecode2];
		$rfsgenename = $genecode2 . "_rfs";
		$rfs[$genecode2] = $_POST[$rfsgenename];

		//get codon positions
		if (!empty($gene_positions2[$genecode2])) {
			$current_gpos = $gene_positions2[$genecode2];
			foreach ( $current_gpos as $k2=> $c2){ //putting choosen codon positions for genes into array in array
				if ($c2 == 'on')	{
					$gene_positions[$genecode2][] =  $k2;
				}
			}
		}
		else { $gene_positions[$genecode2] = array('all'); }
		if (in_array("aas", $gene_positions[$genecode2])) { 
			unset( $gene_positions[$genecode2] ); 
			$gene_positions[$genecode2] = array('aas');
			$gene_by_positions[$genecode2] = 'asone';
		}
		elseif (in_array("all", $gene_positions[$genecode2]) || empty($gene_positions[$genecode2])|| in_array("1st", $gene_positions[$genecode2]) && in_array("2nd", $gene_positions[$genecode2]) && in_array("3rd", $gene_positions[$genecode2])) { 
			unset( $gene_positions[$genecode2] ); 
			$gene_positions[$genecode2] = array('all');
		}
		if (!in_array("aas", $gene_positions[$genecode2]) && !in_array("all", $gene_positions[$genecode2])){
			$gene_by_positions[$genecode2] = 'asone';
		}
		//get partition-by-which-codon-position
		$current_gbypos = $gene_by_positions2[$genecode2];
		// do some test for introns
		if ($ignore_introns == 'no'){
			$query_i = "SELECT intron FROM ". $p_ . "genes WHERE geneCode='$genecode2'";
			$result_i = mysql_query($query_i) or die("Error in query: $query. " . mysql_error());
			// if records present
			if( mysql_num_rows($result_i) > 0 ) {
				while($row_i = mysql_fetch_object($result_i) ) {
					$gene_intr = $row_i->intron;
				}
			}
			// if(!in_array('all', $gene_positions[$genecode2]) || !in_array('aas', $gene_positions[$genecode2]) && $gene_intr != '' && $gene_intr != 'NULL'){
				// $errorList[] = "Must ignore introns to exclude certain codon positions!
											// </br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
											// Set it to 'yes' or chose 'all' for <b>".$genecode2."</b>!";
			// }
			// if ($gene_by_positions[$genecode2] != 'asone' && $gene_intr != '' && $gene_intr != 'NULL'){
				// $errorList[] = "Must ignore introns to partition by certain codon positions!
											// </br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
											// Set it to 'yes' or chose 'asone' for <b>".$genecode2."</b>!";
			// }
		}
	}
		foreach ($geneCodes AS $item) { // building full dataset bp count, setting charset_count[] values, 
									// reading frames and genetic_code
		//$item = clean_item($item);
		$gCquery = "SELECT geneCode, length, readingframe, prot_code, genetic_code, aligned FROM ". $p_ . "genes WHERE geneCode='$item'";
		$gCresult = mysql_query($gCquery) or die("Error in query: $query. " . mysql_error());
		// if records present
		if( mysql_num_rows($gCresult) > 0 ) {
			while( $row = mysql_fetch_object($gCresult) ) {
				$prot_code[$item] = $row->prot_code;
				$genetic_codes[$item] = $row->genetic_code;
				$aligned[$item] = $row->aligned;
				if ($prot_code[$item] != 'yes' && $gene_by_positions[$genecode2] != 'asone'){
					$errorList[] = "The gene <b>$item</b> has been specified</br> as NOT protein coding!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Cannot do else than include all chars as DNA for this!";
				}
				if ($format != 'FASTA'){
					if ($aligned[$item] != 'yes') {
						if ($aligned[$item] == 'no') { 
							$errorList[] = "Gene <b>$item</b> is set as <b>unaligned</b>!
								</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
								Please change that in the gene edit section!</br>
								Will not make a dataset (excl. FASTA) of unaligned data!";
						}
						else { 
							$errorList[] = "Gene <b>$item</b> alignment status is <b>not set</b>!
								</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
								Please change that in the gene edit section!</br>
								Will not make a dataset (excl. FASTA) of unaligned data!";
						}
					}
				}
				if (!in_array('all', $gene_positions[$item])){
					if ($prot_code[$item] == 'yes'){
						if (!in_array($row->genetic_code, array('1','2','3','4','5','6','9','10','11','12','13','14','15','16','21','22','23'))){ 
							$errorList[] = "The <b>genetic code</b> for <b>$item</b> has not been specified!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Please do that in the gene edit section!";
						}
					}
					elseif ($prot_code[$item] == 'notset'){
							$errorList[] = "The gene <b>$item</b> has not been specified</br> as protein coding or not!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Please do that in the gene edit section!";
						}
					elseif ($prot_code[$item] == 'no'){ 
						$errorList[] = "The gene <b>$item</b> has been specified</br> as NOT protein coding!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Cannot do else than include all chars as DNA for this!";
					}
					else {
						$errorList[] = "The gene <b>$item</b> has not been specified correctly!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Please do that in the gene edit section!";
					}
				}
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
				if ( $rf != "1" && $rf != "2" && $rf != "3" && $prot_code[$item] == 'yes'){ 
					if ( $gene_by_positions[$item] == "123" || ! in_array("all", $gene_positions[$item]) || in_array('aas', $gene_positions[$item])) {
						$errors = "Gene $item doesn't have a specified reading frame!
										</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
						if ($gene_by_positions[$item] =="123") { $errors .="Cannot use 12+3 partioning";}
						elseif (in_array('aas', $gene_positions[$item])) { $errors .="Cannot translate to amino acids";}
						else  { $errors .="Cannot use individual position choices";}
						$errorList[] = $errors;
					}
				}
				if ( $gene_by_positions[$item] == "123" || $gene_by_positions[$item] == "each" && ! in_array("all", $gene_positions[$item]) || in_array('aas', $gene_positions[$item])) {
					$errors = "Gene $item doesn't have a specified reading frame!
										</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
										";
										
				}
			}
		}
		else {
			$errorList[] = "Gene <b>$item</b> does not exist in database!</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please add it in the gene edit section!";
		}
	}
	unset($item);
}

// #################################################################################
// Section: Amino acid mode - set variables
// #################################################################################
// checking to see if Amino acid mode is enabled, and in that case just copy the values and fix the by-gene values and proceed to building the dataset
// elseif (isset($_POST['geneCodes2']) && isset($_POST['codes2']) && isset($_POST['genetic_codes'])){
	// $format = $_POST['format2']; //echo "format = >$format<</br>";
	// $geneCodes = explode(",", $_POST['geneCodes2']);
	// $codes = explode(",", $_POST['codes2']); // $codes3 = implode(",",$codes);echo "codes: " . $codes3 ."</br>";
	// $taxonadds = explode(",",$_POST['taxonadds2']); //$taxonadds3 = implode(",",$taxonadds);echo "taxonadds: " . $taxonadds3 ."</br>";
	// $outgroup = $_POST['outgroup2']; //echo "outgroup: $outgroup</br>";
	// $positions = array("aas");
	// $by_positions = "asone";
	// $number_of_taxa = $_POST['number_of_taxa'];
	// $genetic_codes = $_POST['genetic_codes'];
	// $ignore_introns = $_POST['ignore_introns'];
	// foreach ($geneCodes as $genecode2){
		// //get charset counts and rfs
		// $charset_count[$genecode2] = $_POST[$genecode2];
		// $rfsgenename = $genecode2 . "_rfs";
		// $rfs[$genecode2] = $_POST[$rfsgenename];
	// }
// }

// #################################################################################
// Section: normal mode - set variables
// #################################################################################
 else { //if special mode or Amino acid mode is not enabled, checking and building from the beginning
	if (isset($_POST['geneCodes'])){
		foreach ( $_POST['geneCodes'] as $k1=> $c1){ //putting choosen genes into array
			if ($c1 == 'on')	{
				$geneCodes[] =  $k1;
			}
		}
	}
	else {
		$errorList[] = "No genes choosen - Please try again!"; 
	}

	 //input data
	if (trim($_POST['codes']) != "") {
		$raw_codes = explode("\n", $_POST['codes']);
	}
	else { 
		unset($raw_codes); 
	}

	$format = $_POST['format'];
	if ( !isset($format) || $format == '') {
		$errorList[] = "No dataset FORMAT choosen - Please try again!";
	}

	$by_positions = $_POST['by_positions'];

	if ($format == "NEXUS" || $format == "TNT" || $format == "PHYLIP" ) {
		if (isset($_POST['outgroup']) && trim($_POST['outgroup']) != "") {
			$outgroup = clean_item($_POST['outgroup']);
			$outgroup = trim($outgroup);
			//$outgroup = trim($_POST['outgroup']);
		}
		else {
			unset($outgroup);
		}
	}
	else {
		unset($outgroup);
	}

	$taxonadds = array();
	foreach ( $_POST['taxonadds'] as $k=> $c){
		if ($c == 'on')	{
			$taxonadds[] = $k;
		}
	}
	$exclude_missing = $_POST['exclude_missing'];
	$less_than_genes = trim($_POST['less_than_genes']);
	$ignore_introns = trim($_POST['ignore_introns']);

	// if ($format != "FASTA") { removing this - thus allowing for fasta retrieval of certain positions
	if ( isset($_POST['positions'])){
		$positions = array();
		foreach ( $_POST['positions'] as $k2=> $c2){ //putting choosen genes into array
			if ($c2 == 'on')	{
				$positions[] =  $k2;
			}
		}
	}else {$positions = array("all"); }

	// do some test for "position choices"
	if (in_array('aas', $positions)){ // setting up for amino acid partition mode
		unset( $positions ); 
		$positions = array('aas') ;
		$by_positions = "asone";
		}
	elseif (in_array('special', $positions)){ // setting up for special gene codon partition mode
		unset( $positions ); 
		$positions = array('special') ;
	}
	elseif (in_array("all", $positions) || in_array("1st", $positions) && in_array("2nd", $positions) && in_array("3rd", $positions)) { 
		unset( $positions ); 
		$positions = array('all') ;
	}

	elseif ( ! in_array("all", $positions) && count($positions) < 2 && $by_positions !='special') { //if only one codon position choosen - force "asone"
		$by_positions = "asone";
	}
	else {
		if ($by_positions == "123") {
			$errorList[] = "Cannot use the '12+3' partitioning without including all positions!";
		}
	}
	// do some test for introns
	// if ($ignore_introns == 'no'){
		// if( !in_array('special', $positions) && !in_array('all', $positions)){
			// $errorList[] = "Must ignore introns to exclude certain codon positions!
										// </br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
										// Set it to 'yes'!";
		// }
		// if ($by_positions != 'asone'){
			// $errorList[] = "Must ignore introns to partition by certain codon positions!
										// </br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
										// Set it to 'yes'!";
		// }
	// }

	$result = mysql_query("set names utf8") or die("Error in query: $query. " . mysql_error());

	$seq_string = "?"; #this will be a replacement for NULL sequences
	$bp = 0; #number of base pairs

	//check outgroup existance
	if ($format == "NEXUS" || $format == "TNT" || $format == "PHYLIP"){
		if (isset($outgroup)) {
			$queryog = "SELECT code, genus, species, orden, family, subfamily, tribe, subtribe, subspecies, hostorg, auctor FROM ". $p_ . "vouchers WHERE code='$outgroup'";
			$resultog = mysql_query($queryog) or die("Error in query: $query. " . mysql_error());
			// if records present
			if( mysql_num_rows($resultog) > 0 ) {
				$outgroup_to_taxa = $outgroup;
				if ($format == "NEXUS"){
					while( $rowog = mysql_fetch_object($resultog) ) {
						$currentcode = $rowog->code;
							foreach ( $_POST['taxonadds'] as $k=> $c){
								if ($c == 'on')	{
									$taxarray[] .=  $rowog->$k;
								}
							}
							$taxon .= implode("_" , $taxarray);
							$replaces = array(" ","________","_______","______","_____","____","___","__");
							$taxon = str_replace($replaces, "_", $taxon);
							$replaces2 = array("(" , ")" , ";" , ",", "=", "?", "\"", "/");
							$taxon = str_replace($replaces2, "", $taxon);
							$taxon = str_replace("-", "_", $taxon);
							$taxon = str_replace($replaces, "_", $taxon);
							$taxon = substr($taxon, 0, 75);
							$taxon = "$taxon ";
							$outgroup = rtrim($taxon);
					}
				}
			}
			else { $errorList[] = "The specified outgroup: <b>$outgroup</b> does not exist in db!";}
		}
	}
	foreach ($geneCodes AS $item) { // building full dataset bp count, setting charset_count[] values, 
									// reading frames and genetic_code
		//$item = clean_item($item);
		$gCquery = "SELECT geneCode, length, readingframe, prot_code, genetic_code, aligned FROM ". $p_ . "genes WHERE geneCode='$item'";
		$gCresult = mysql_query($gCquery) or die("Error in query: $query. " . mysql_error());
		// if records present
		if( mysql_num_rows($gCresult) > 0 ) {
			while( $row = mysql_fetch_object($gCresult) ) {
				$prot_code[$item] = $row->prot_code;
				$genetic_codes[$item] = $row->genetic_code;
				$aligned[$item] = $row->aligned;
				if ($format != 'FASTA'){
					if ($aligned[$item] != 'yes') {
						if ($aligned[$item] == 'no') { 
							$errorList[] = "Gene <b>$item</b> is set as <b>unaligned</b>!
								</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
								Please change that in the gene edit section!</br>
								Will not make a dataset (excl. FASTA) of unaligned data!";
						}
						else { 
							$errorList[] = "Gene <b>$item</b> alignment status is <b>not set</b>!
								</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
								Please change that in the gene edit section!</br>
								Will not make a dataset (excl. FASTA) of unaligned data!";
						}
					}
				}
				if (!in_array('special', $positions) && !in_array('all', $positions)){
					if ($prot_code[$item] == 'yes'){
						if (!in_array($row->genetic_code, array('1','2','3','4','5','6','9','10','11','12','13','14','15','16','21','22','23'))){ 
							$errorList[] = "The <b>genetic code</b> for <b>$item</b> has not been specified!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Please do that in the gene edit section!";
						}
					}
					elseif ($prot_code[$item] == 'notset'){
							$errorList[] = "The gene <b>$item</b> has not been specified</br> as protein coding or not!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Please do that in the gene edit section!";
						}
					elseif ($prot_code[$item] == 'no'){ 
						unset($genetic_codes[$item]);
					}
					else {
						$errorList[] = "The gene <b>$item</b> has not been specified correctly!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
									Please do that in the gene edit section!";
					}
				}
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
				if ( $rf != "1" && $rf != "2" && $rf != "3" && $prot_code[$item] == 'yes'){ 
					if ( $by_positions == "123" || ! in_array("all", $positions) || in_array('aas', $positions)) {
						$errors = "Gene $item doesn't have a specified reading frame!
										</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
						if ($by_positions =="123") { $errors .="Cannot use 12+3 partioning";}
						elseif (in_array('aas', $positions)) { $errors .="Cannot translate to amino acids";}
						else  { $errors .="Cannot use individual position choices";}
						$errorList[] = $errors;
					}
				}
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

	//setting outgroups as first taxa in list for TNT and PHY datasets
	if (isset($outgroup_to_taxa)) {
		unset($codes_wo_og);
		if ( $format == "TNT" || $format == "PHYLIP") {
			if (in_array($outgroup_to_taxa, $codes ) ) { 
				$codes_wo_og = array();
				foreach ($codes as $code1) { 
					if ($code1 != $outgroup_to_taxa){
						$codes_wo_og[] = $code1;
					}
				} 
			}
		}
		if (isset($codes_wo_og)){ $codes = $codes_wo_og;}
		if ( ! in_array($outgroup_to_taxa, $codes) ) { 
			// $errorList[] = "The specified outgruop: <b>$outgroup</b> does not exist among the dataset codes!";
			array_unshift( $codes, $outgroup_to_taxa );
		}
	}
	
	// for single genes if exclude taxa box "yes" is checked, removing taxa that lacks that gene sequence
	if (count($geneCodes) == 1 && $exclude_missing == "yes"){
	$codes_to_remove = array();
		foreach($codes AS $item) {
			$cquery = "SELECT sequences FROM ". $p_ . "sequences WHERE code='$item' AND geneCode='$geneCodes[0]'";
			$cresult = mysql_query($cquery) or die("Error in query: $query. " . mysql_error());
			// if records present
			if( mysql_num_rows($cresult) > 0 ) {
				while( $row = mysql_fetch_object($cresult) ) {
					$ctrlseq = clean_item(trim($row->sequences));
					$del_list = array("'\?'", "'\-'", "'N'", "'X'", "'\s'");
					foreach ($del_list as $del_char) {$ctrlseq = preg_replace($del_char,'',$ctrlseq);}
					if ($ctrlseq == '') {$codes_to_remove[] = $item;}
				}
			}
			else {$codes_to_remove[] = $item;}
		}
		if (count($codes_to_remove) > 0) {$codes = array_diff($codes, $codes_to_remove);}
	}
	// for multi gene datasets with specified minimum number of genes present, remove taxa with less
	if ($less_than_genes != '' && count($geneCodes) > 1){
		if (!is_numeric($less_than_genes) && $less_than_genes != '') {$errorList[] = "Multigene taxa removal gene count number is not in numeric format!</br>
										&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
										Please add a number between 1 and max number of genes!";
		}
		else {
			$less_than_genes = intval($less_than_genes);
			if ( $less_than_genes > 0 && $less_than_genes <= count($geneCodes)) {
				$less_than_genes = $less_than_genes;
				$codes_to_remove = array();
				$genecodeloop = array();
				foreach($codes AS $item) {
					$ng = count($geneCodes);
					foreach ($geneCodes as $gcl) {
						$cquery = "SELECT sequences FROM ". $p_ . "sequences WHERE code='$item' AND geneCode='$gcl'";
						$cresult = mysql_query($cquery) or die("Error in query: $query. " . mysql_error());
						// if records present
						if( mysql_num_rows($cresult) > 0 ){
							while( $row = mysql_fetch_object($cresult) ) {
								$ctrlseq = clean_item(trim($row->sequences));
								$del_list = array("'\?'", "'\-'", "'N'", "'X'", "'\s'");
								foreach ($del_list as $del_char) {$ctrlseq = preg_replace($del_char,'',$ctrlseq);}
								if ($ctrlseq == '') {$ng = $ng - 1; }
							}
						}
						else {$ng = $ng - 1;}
					}
				if ($ng < $less_than_genes) {$codes_to_remove[] = $item;}
				}
				if (count($codes_to_remove) > 0) {$codes = array_diff($codes, $codes_to_remove);}
			}
			elseif ($less_than_genes > count($geneCodes)) {
				$errorList[] = "Trying to remove taxa with more genes</br>
								than total number of genes!!
								</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
								Please enter a number from 1 to max number of genes!";
			}
		}
	}

	
	//if no taxa left or presented
	$number_of_taxa = count($codes);
	if ($number_of_taxa == 0) {$errorList[] = "No codes specified! No use creating empty datasets...
												</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
												Please go back and add voucher codes to run!"; 
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
// Section: Start building dataset or specials or amino acid choice outputs
// #################################################################################

// #################################################################################
// Section: Creating special mode choice output
// #################################################################################
	if ( in_array("special", $positions) && !isset($gene_positions2) && !isset($gene_by_positions2) ){
			$title = "$config_sitename: Dataset Special";
		// print html headers
		$admin = false;
		$in_includes = true;
		include_once 'header.php';
		// print navegation bar
		nav();
		// begin HTML page content
		echo "<div id=\"content_narrow\">";
		echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
				<tr><td valign=\"top\" colspan=\"2\"><h1>Create special dataset</h1>
				<p>Enter the required info to make yourself a ready-to-run dataset in $format format:<br />";
		// print as list
		echo '<form action="process_dataset.php" method="post">';
		echo "Choose wanted codon positions for the separate genes ('all' will override other choices):";
		echo '<br><ul></td></tr>';
		foreach ($geneCodes as $genes) {
			echo "<td>Gene $genes: </td><td>";
			echo "<input type=\"checkbox\" name=\"gene_positions2[$genes][all]\" checked>all&nbsp;&nbsp;&nbsp;";
			echo "<input type=\"checkbox\" name=\"gene_positions2[$genes][1st]\" >1st&nbsp;&nbsp;&nbsp;";
			echo "<input type=\"checkbox\" name=\"gene_positions2[$genes][2nd]\" >2nd&nbsp;&nbsp;&nbsp;";
			echo "<input type=\"checkbox\" name=\"gene_positions2[$genes][3rd]\" >3rd&nbsp;&nbsp;&nbsp;";
			echo "<input type=\"checkbox\" name=\"gene_positions2[$genes][aas]\" >Amino acids&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;";
			echo "<input type=\"radio\" name=\"gene_by_positions2[$genes]\" value=\"asone\" checked>as one&nbsp;&nbsp;&nbsp;";
			echo "<input type=\"radio\" name=\"gene_by_positions2[$genes]\" value=\"each\">each&nbsp;&nbsp;&nbsp;";
			echo "<input type=\"radio\" name=\"gene_by_positions2[$genes]\" value=\"123\">12+3&nbsp;&nbsp;&nbsp;";
			echo "</td></tr>";
		}
		//keeping old values
		$geneCodes2 = implode(",", $geneCodes);
		$codes2 = implode(",", $codes);
		$positions2 = implode(",", $positions);
		$taxonadds2 = implode(",", $taxonadds);
		?>
			<input type="hidden" name="format2" value="<?php echo $format; ?>" >
			<input type="hidden" name="geneCodes2" value="<?php echo $geneCodes2; ?>" >
			<input type="hidden" name="codes2" value="<?php echo $codes2; ?>" >
			<input type="hidden" name="taxonadds2" value="<?php echo $taxonadds2; ?>" >
			<input type="hidden" name="ignore_introns" value="<?php echo $ignore_introns; ?>" >
			<?php if (isset($outgroup)) { ?> <input type="hidden" name="outgroup2" value="<?php echo $outgroup; ?>" > <?php } ?>
			<input type="hidden" name="number_of_taxa" value="<?php echo $number_of_taxa; ?>" >
			<?php if (isset($prot_code)) { ?> <input type="hidden" name="prot_code" value="<?php echo $prot_code; ?>" > <?php } ?>
			<?php if (isset($genetic_codes)) { ?> <input type="hidden" name="genetic_codes" value="<?php echo $genetic_codes; ?>" > <?php } ?>
			<?php if (isset($aligned)) { ?> <input type="hidden" name="aligned" value="<?php echo $aligned; ?>" > <?php } ?>
			<?php foreach ($geneCodes as $genes){
				$value = $charset_count[$genes];
				$rfs2 = $rfs[$genes];
				$rfsgenename = $genes . "_rfs";
				echo "<input type='hidden' name='$genes' value='$value' >";
				echo "<input type='hidden' name='$rfsgenename' value='$rfs2' >";
			}unset($genes, $value); 
			?>
			</td></tr><tr><td><input type="submit" name="process_dataset" value="Continue dataset creation" /></td></tr>
			</form>
			</tr>
				</table> <!-- end super table -->
				</div> <!-- end content -->
		<?php
		//make footer
		make_footer($date_timezone, $config_sitename, $version, $base_url);
		echo '</body></html>';
	}

	// end specials choose list --------------------------------------------------------------------------------------------

// #################################################################################
// Section: Creating amino acid choice output
// #################################################################################
	// elseif ( in_array("aas", $positions) && !isset($genetic_codes) ){
			// $title = "$config_sitename: Amino Acid Dataset";
		// print html headers
		// $admin = false;
		// $in_includes = true;
		// include_once 'header.php';
		// print navegation bar
		// nav();
		// begin HTML page content
		// echo "<div id=\"content_narrow\">";
		// echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
				// <tr><td valign=\"top\" colspan=\"2\"><h1>Create Amino acid dataset</h1>
				// <p>Enter the required info to make yourself a ready-to-run dataset in $format format:<br />";
		// print as list
		// echo '<form action="process_dataset.php" method="post">';
		// echo "Choose genetic code for translation:";
		// echo '<br><ul></td></tr>';
		// foreach ($geneCodes as $genes) {
			// echo "<td>Gene $genes: </td><td>";
			// ? >
			// <select name=<?php echo "genetic_codes[$genes]";? > size="1" style=" BORDER-BOTTOM: outset; BORDER-LEFT: 
			// outset; BORDER-RIGHT: outset; BORDER-TOP: outset; FONT-FAMILY: 
			// Arial; FONT-SIZE: 12px"> 
			  // <!-- create a pulldown-list with all taxon set names in the db -->
			    // <option value=1 selected>Standard
                // <option value=2>Vertebrate Mitochondrial
                // <option value=3>Yeast Mitochondrial
                // <option value=4>Mold, Protozoan and Coelenterate Mitochondrial. Mycoplasma, Spiroplasma
                // <option value=5>Invertebrate Mitochondrial
                // <option value=6>Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
                // <option value=9>Echinoderm Mitochondrial
                // <option value=10>Euplotid Nuclear
                // <option value=11>Bacterial and Plant Plastid
                // <option value=12>Alternative Yeast Nuclear
                // <option value=13>Ascidian Mitochondrial
                // <option value=14>Flatworm Mitochondrial
                // <option value=15>Blepharisma Macronuclear
                // <option value=16>Chlorophycean Mitochondrial
                // <option value=21>Trematode Mitochondrial
                // <option value=22>Scenedesmus obliquus mitochondrial
                // <option value=23>Thraustochytrium mitochondrial code
				// </select></td></tr>
				// <?php
		// }
		// keeping old values
		// $geneCodes2 = implode(",", $geneCodes);
		// $codes2 = implode(",", $codes);
		// $positions2 = implode(",", $positions);
		// $taxonadds2 = implode(",", $taxonadds);
		// ? >
			// <input type="hidden" name="format2" value="<?php echo $format; ? >" >
			// <input type="hidden" name="geneCodes2" value="<?php echo $geneCodes2; ? >" >
			// <input type="hidden" name="codes2" value="<?php echo $codes2; ? >" >
			// <input type="hidden" name="taxonadds2" value="<?php echo $taxonadds2; ? >" >
			// <?php if (isset($outgroup)) { ? > <input type="hidden" name="outgroup2" value="<?php echo $outgroup; ? >" > <?php } ? >
			// <input type="hidden" name="number_of_taxa" value="<?php echo $number_of_taxa; ? >" >
			// <input type="hidden" name="ignore_introns" value="<?php echo $ignore_introns; ? >" >
			// <?php foreach ($geneCodes as $genes){
				// $value = $charset_count[$genes];
				// $rfs2 = $rfs[$genes];
				// $rfsgenename = $genes . "_rfs";
				// echo "<input type='hidden' name='$genes' value='$value' >";
				// echo "<input type='hidden' name='$rfsgenename' value='$rfs2' >";
			// }unset($genes, $value); 
			// ? >
			// </td></tr><tr><td><input type="submit" name="process_dataset" value="Continue dataset creation" /></td></tr>
			// </form>
			// </tr>
				// </table> <!-- end super table -->
				// </div> <!-- end content -->
		// < ?php
		// make footer
		// make_footer($date_timezone, $config_sitename, $version, $base_url);
		// echo '</body></html>';
	//}
	// end aa translation code choose list --------------------------------------------------------------------------------------------
	// #################################################################################
	// Section: Build dataset and adding choosen name-extensions
	// #################################################################################
	else {
		$num_genes = 0;
		$output_lines = "";
		$taxout_array = array();
		$seqout_array = array();
		$intron_dataset = array();
		$intron_lengths = array();
		foreach ($geneCodes AS $geneCode) {
			$num_genes = $num_genes + 1;
			$aa_or_not[$geneCode] = 'no';
			// check for introns
			//if ($ignore_introns == 'yes'){
				$query_i = "SELECT intron, genetic_code, prot_code FROM ". $p_ . "genes WHERE geneCode='$geneCode'";
				$result_i = mysql_query($query_i) or die("Error in query: $query. " . mysql_error());
				// if records present
				if( mysql_num_rows($result_i) > 0 ) {
					while($row_i = mysql_fetch_object($result_i) ) {
						$gene_int = $row_i->intron;
						if (isset($gene_int) && $gene_int != ''&& $gene_int != 'NULL'){
							$intron = remove_introns(str_pad("A", $charset_count[$geneCode], "A"),$gene_int );
							$new_length = $intron[1];
							$number_of_introns = $intron[2];
							for ($i = 1; $i <= $intron[2]; $i++){
								$intron_lengths[$geneCode][$i] = strlen($intron[3+$i]);
							}
						}else $new_length = 'no';
					}
				}
				else {$new_length = 'no';}
			//}
			foreach ($codes AS $item) {
				$item = clean_item($item);
				$query = "SELECT code, genus, species, orden, family, subfamily, tribe, ";
				$query .= "subtribe, subspecies, hostorg, auctor FROM ". $p_ . "vouchers WHERE code='$item'";
				$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
				// if records present
				if( mysql_num_rows($result) > 0 ) {
					while( $row = mysql_fetch_object($result) ) {
						$seq = "?";
						$currentcode = $row->code;
						// #################################################################################
						// Section: Name builder
						// #################################################################################
						$taxon = "";
						$taxarray = array();
						if( $format == "FASTA" ) {	
							$taxon .= ">";	
						}
						//foreach ( $_POST['taxonadds'] as $k=> $c){
						//	if ($c == 'on')	{
						foreach ($taxonadds as $k){
							if ($k == genecode) { 
								if ($format == "NEXUS") {
									$taxarray[] .=  "[$geneCode]";
								} 
								else {
									$taxarray[] .=  "$geneCode";	
								}
							}
							else {	
								$taxarray[] .=  $row->$k;	
							}
						}
						//	}
						//}
						$taxon .= implode("_" , $taxarray);
						$replaces = array(" ","________","_______","______","_____","____","___","__" );
						$taxon = str_replace($replaces, "_", $taxon);
						$taxon = str_replace(">_", ">", $taxon);
						if( $format != "FASTA" ) {
							$replaces2 = array("(" , ")" , ";" , ",", "=", "?", "\"", "/");
							$taxon = str_replace($replaces2, "", $taxon);
							$taxon = str_replace("-", "_", $taxon);
							$taxon = str_replace($replaces, "_", $taxon);
							if ( $format == "NEXUS" ) { 
								//$taxon = substr($taxon, 0, 75); $taxon = "'$taxon'";
								$taxon = substr(str_pad($taxon, 55, " "), 0, 55);
							}
							else {
								$taxon = substr(str_pad($taxon, 55, " "), 0, 55);
							}
						}
						// if (in_array("geneCode", $_POST['taxonadds'])) {$taxon = str_replace("_[", "[", $taxon);}
						if (in_array("geneCode", $taxonadds)) {$taxon = str_replace("_[", "[", $taxon);}
						if ($format == "PHYLIP" && $num_genes > 1) {	$taxon = str_pad("", 55, " ");	}
						if ($format != "FASTA" ) { $taxon = "$taxon ";}
						$taxout_array[$geneCode][$item] = $taxon;
							
						// #################################################################################
						// Section: Sequence builder
						// #################################################################################
						$query_b = "SELECT sequences FROM ". $p_ . "sequences WHERE code='$row->code' AND geneCode='$geneCode'";
						$result_b = mysql_query($query_b) or die("Error in query: $query_b. " . mysql_error());
						// if records present
						if( mysql_num_rows($result_b) > 0 ) {
							while( $row_b = mysql_fetch_object($result_b) ) {
								$seq = $row_b->sequences;
								//if( $format == "FASTA" ) { // do nothing - just present raw sequences
										//$seq = "\n" . $seq;
								if ($ignore_introns == 'yes' && $new_length != 'no'){
									$seq_i = remove_introns($seq, $gene_int);
									$seq = $seq_i[0];
								}
								if ($ignore_introns == 'no' && $new_length != 'no' && $format != "FASTA" && $format != "TNT"){
									$seq_i = remove_introns($seq, $gene_int);
									$seq = $seq_i[0];
									for ($i = 1; $i <= $seq_i[2]; $i++){
										$intron_dataset[$geneCode]["_intron_". $i][$item] = $seq_i[$i+3];
									}
								}
								$seqout_array[$geneCode][$item] = $seq;
								// }
								//else {						// do something
								if ( isset($gene_positions)){ $positions = $gene_positions[$geneCode];} // if special mode
								if ($format!="FASTA" && strlen($charset_count[$geneCode] < strlen($seq))){ // checking for too long sequence
									$errorList[] = "The $geneCode sequence of $item is longer (". strlen($seq) . ">" . $charset_count[$geneCode] .")that the specified gene length!
									</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please edit gene length or check the sequence";
								}
								elseif (in_array("aas", $positions) &&  isset($genetic_codes) && $prot_code[$geneCode] == 'yes') {
									// translating amino acid sequence and replacing ? and saces with X
									if ($rfs[$geneCode] == "2") { $seq = substr($seq,1);}
									elseif ($rfs[$geneCode] == "3") { $seq = substr($seq,2);}
									$seq = translate_DNA_to_protein($seq,$genetic_codes[$geneCode]);
									//$seq = str_replace("?", "X", $seq);
									$seq = str_replace("X", "?", $seq);
									$seq = str_replace(" ", "?", $seq);
									$seqout_array[$geneCode][$item] = $seq;
									$aa_or_not[$geneCode] = 'yes';
									//$seqout_1 = translate_DNA_to_protein($seq,$genetic_codes[$geneCode]);
									//echo ">$item</br>$seqout_1</br>";
								}
								else {
									//if ( isset($gene_positions)){ $positions = $gene_positions[$geneCode];} // if special mode
									// create choosen codon position sequences
									if (! in_array("all", $positions) && $prot_code[$geneCode] == 'yes') {
										if ($rfs[$geneCode] == "1") { $num_nuc = "1";}
										elseif ($rfs[$geneCode] == "2") { $num_nuc = "3";}
										elseif ($rfs[$geneCode] == "3") { $num_nuc = "2";}
										else { $num_nuc = "1";}
										$pos_array = array();
										$sequence_array = preg_split('#(?<=.)(?=.)#s', $seq); // making sequence/nucleotide array
										foreach ($sequence_array as $nuc){
											$pos_array[$num_nuc][] = $nuc;
											if ($num_nuc == 1 || $num_nuc == 2) { $pos_array[12][] = $nuc;}
											if ($num_nuc == 1 || $num_nuc == 3) { $pos_array[13][] = $nuc;}
											if ($num_nuc == 2 || $num_nuc == 3) { $pos_array[23][] = $nuc;}
											if ($num_nuc == 3 ) { $num_nuc = 1;}
											else {$num_nuc = $num_nuc + 1;}
											array_shift($sequence_array);
										}
										if ( in_array("1st", $positions) && ! in_array("2nd", $positions) && ! in_array("3rd", $positions)) { $seqout_array[$geneCode][$item] = implode($pos_array[1]);}
										elseif (! in_array("1st", $positions) && in_array("2nd", $positions) && ! in_array("3rd", $positions)) { $seqout_array[$geneCode][$item] = implode($pos_array[2]);}
										elseif (! in_array("1st", $positions) && ! in_array("2nd", $positions) && in_array("3rd", $positions)) { $seqout_array[$geneCode][$item] = implode($pos_array[3]);}
										elseif ( in_array("1st", $positions) && in_array("2nd", $positions) && ! in_array("3rd", $positions)) { $seqout_array[$geneCode][$item] = implode($pos_array[12]);}
										elseif ( in_array("1st", $positions) && ! in_array("2nd", $positions) &&  in_array("3rd", $positions)) { $seqout_array[$geneCode][$item] = implode($pos_array[13]);}
										elseif (! in_array("1st", $positions) && in_array("2nd", $positions) &&  in_array("3rd", $positions)) { $seqout_array[$geneCode][$item] = implode($pos_array[23]);}
										else { $seqout_array[$geneCode][$item] = $seq; }
									}
									else { $seqout_array[$geneCode][$item] = $seq; }
									// padding string for aligned datasets
									//$seq = str_pad($seq, $charset_count[$geneCode], "?");
								}
								// }
							}
						}
						else {
							if( $format == "FASTA" && count($geneCodes) == 1) {
								$seq = "\n";
							}
							else {
								if ($ignore_introns == 'no' && $new_length != 'no' && $format != "FASTA" && $format != "TNT"){
									$seq_i = remove_introns($seq, $gene_int);
									$seq = $seq_i[0];
									for ($i = 1; $i <= $seq_i[2]; $i++){
										$intron_dataset[$geneCode]["_intron_". $i][$item] = str_pad("?", $intron_lengths[$geneCode][$i], "?");
									}
								}
								if (in_array("aas", $positions)) { 
									$seq = "?"; 
								}
								else { 
									$seq = "?"; 
									//$seq = str_pad($seq, $charset_count[$geneCode], "?");
								}
								$seqout_array[$geneCode][$item] = $seq;
							}
						}
							//$output_lines .= $seq . "\n";
					}
				}
			}
		}
		unset($item);
	
		// #################################################################################
		// Section: setting bp numbers for partitions if needed
		// #################################################################################
		if ($format == "FASTA" && count($geneCodes) == 1 ) { }//&& ! in_array("all", $positions)) {
		else{
			unset ($charset_count);
			$datatype_mixed = array();
			$datatype_mixed_sum = 1;
			if (isset($seqout_array)){
				foreach ($seqout_array as $g => $s) { 
					if (in_array("aas", $positions) || $aa_or_not[$g] == 'yes') {$charset_count[$g] = 0; }
					foreach ($s as $n) {
						if (strlen($n) > $charset_count[$g] ) { $charset_count[$g] = strlen($n);}
					}
					if ($aa_or_not[$g] == 'yes') { $datatype_mixed[] = "protein:$datatype_mixed_sum-". ($charset_count[$g]+$datatype_mixed_sum-1);}
					else {$datatype_mixed[] = "dna:$datatype_mixed_sum-". ($charset_count[$g]+$datatype_mixed_sum-1);}
					$datatype_mixed_sum = $charset_count[$g]+$datatype_mixed_sum;
					if (count($intron_dataset[$g])>0 && $format != "FASTA" && $format != "TNT" && $ignore_introns == 'no') {
							$i = 1;
							foreach ($intron_lengths[$g] as $il) {
								$charset_count[$g."_".$i] = $il;
								$datatype_mixed[] = "dna:$datatype_mixed_sum-". ($il+$datatype_mixed_sum-1);
								$datatype_mixed_sum = $il+$datatype_mixed_sum;
								$i++;
							}
					}
				} unset($s, $g, $n);
			$bp = array_sum($charset_count);
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
		else{ //start building dataset
			
			// #################################################################################
			// Section: Build output - intro lines
			// #################################################################################
			if( $format == "TNT" ) {  // creating intro lines
				if (in_array("aas", $positions)) {
					$output = "nstates prot;\nxread\n$bp $number_of_taxa\n";
				}
				else {
					$output = "nstates dna;\nxread\n$bp $number_of_taxa\n";
				}
			}
			elseif( $format == "PHYLIP" ) {
				$output = "$number_of_taxa $bp\n";
				$phy_partitions = array();
			}
			elseif( $format == "FASTA" ) {
				$output = "";
			}
			else {
				$which_nex_partitions = array();
				$nex_partitions = array();
				if (in_array("aas", $positions) || in_array("yes", $aa_or_not)) {
					if (!in_array("no", $aa_or_not) && $ignore_introns == 'yes'){
						$output = "#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=$number_of_taxa NCHAR=$bp;\nFORMAT INTERLEAVE DATATYPE=PROTEIN MISSING=? GAP=-;\nMATRIX\n";
					}
					else {
						$output = "#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=$number_of_taxa NCHAR=$bp;\nFORMAT INTERLEAVE DATATYPE=";
						$output .= "mixed(". implode(",",$datatype_mixed) .") ";
						$output .= "MISSING=? GAP=-;\nMATRIX\n";
					}
				}
				else {
					$output = "#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=$number_of_taxa NCHAR=$bp;\nFORMAT INTERLEAVE DATATYPE=DNA MISSING=? GAP=-;\nMATRIX\n";
				}
			}
			// #################################################################################
			// Section: Build output - sequence blocks
			// #################################################################################
				
			// creating sequence blocks / gene
			$multigenefasta = array();
			if ($format == "FASTA" && count($geneCodes) > 1){
				foreach ($codes AS $item) {
					$item = clean_item($item);
					$output .= $taxout_array[$geneCodes[0]][$item] . "\n";
					foreach ($geneCodes as $genCod) {
						$output .= str_pad($seqout_array[$genCod][$item], $charset_count[$genCod], "?");
					}
					$output .= "\n";
				}
			}
			else {
				foreach ($geneCodes as $geCo) { 
					## check if there is any sequence at all for this gene, for this list of specimens
					## has_seqs?
					if( has_seqs($seqout_array, $geCo) == "false" ) {
						$warning[] = "Warning, no sequences in partition $geCo"; 
					}

					if( $format == "TNT" ) {
						if (in_array("aas", $positions)) {$output .= "\n&[PROTEIN]\n";}
						else{ $output .= "\n&[dna]\n"; }
					}
					elseif( $format == "FASTA" && count($geneCodes) == 1) {	
						$output .= ">$geCo\n--------------\n";
					}
					elseif ( $format == "NEXUS" ) {
						## check if there is any sequence at all for this gene, for this list of specimens
						## has_seqs?
						$has_seqs = has_seqs($seqout_array, $geCo);
						if( $has_seqs == "true" ) {
							$output .= "\n[$geCo]\n";
						}
						else {
							$output .= "\n[$geCo - warning no sequences in this partition]\n";
						}
					}
					else {}
					foreach ($codes AS $item) {
						$item = clean_item($item);
						$output .= $taxout_array[$geCo][$item];
						if ($format != "FASTA"){
							$output .= str_pad($seqout_array[$geCo][$item], $charset_count[$geCo], "?") . "\n"; 
						}
						else {
							$output .= "\n" . $seqout_array[$geCo][$item] . "\n"; 
							}
						}
					if ($format == 'PHYLIP') {$output .= "\n";}
					if (count($intron_dataset[$geCo]) > 0 && $ignore_introns == 'no' && $format != "FASTA" && $format != "TNT") {
						$Intr = 1;
						foreach ($intron_dataset[$geCo] as $geIn){
							if ($format == "NEXUS") {$output .= "\n[". $geCo . "_intron$Intr" . "]\n";}
							foreach ($codes AS $item) {
								$output .= $taxout_array[$geCo][$item];
								$output .= $geIn[$item] . "\n";
							}
							$Intr++;
						}
						if ($format == 'PHYLIP') {$output .= "\n";}
					}
				}
			}
	// #################################################################################
	// Section: Build output - partition specifying block
	// #################################################################################
	//creating NEXUS and PHYLIP partition-file blocks (for NEXUS only blocks for codon position partitioned
	if( $format == "PHYLIP" || $format == "NEXUS") {
		$phybp = "1";
		$capture_positions = $positions;
		$capture_by_positions = $by_positions;
		foreach ($geneCodes as $gCPHY){
			$positions = $capture_positions;
			$by_positions = $capture_by_positions;
			if ($prot_code[$gCPHY] == 'no') {$positions = array('all'); $by_positions = 'asone';}
			if (isset($gene_positions) && isset($gene_by_positions)){ // if special mode
				$positions = $gene_positions[$gCPHY]; 
				$by_positions = $gene_by_positions[$gCPHY];
			}
			if (count($positions) == '1' && !in_array("all", $positions)) {
				$by_positions = "asone";
			}
			//setting frequency of codon positions
			if ( in_array("all", $positions)) { $p_jumps = "\\3" ;}
			else {  if (count($positions) > 1) {$p_jumps = "\\" . count($positions);} else { $p_jumps = ""; } }
			//setting partitionname ending for "all" sequences
			$part_name_end = "";
			if (in_array("1st", $positions)) { $part_name_end .= "1";}
			if (in_array("2nd", $positions)) { $part_name_end .= "2";}
			if (in_array("3rd", $positions)) { $part_name_end .= "3";}
			//setting readingframes
			if ($rfs[$gCPHY] == "2") { 
				if (in_array("all", $positions)) {
					$p_pos1 = "1";$p_pos2 = "2";$p_pos3 = "0";
				} 
				else { 
					if ( ! in_array("1st", $positions)) { $p_pos2 = "1";$p_pos3 = "0";}
					elseif ( ! in_array("2nd", $positions)) { $p_pos1 = "1";$p_pos3 = "0";}
					elseif ( ! in_array("3rd", $positions)) { $p_pos1 = "0";$p_pos2 = "1";}
				}
			}
			elseif ($rfs[$gCPHY] == "3") {
				if (in_array("all", $positions)) {
					$p_pos1 = "2";$p_pos2 = "0";$p_pos3 = "1";
				}
				else { 
					if ( ! in_array("1st", $positions)) { $p_pos2 = "0";$p_pos3 = "1";}
					elseif ( ! in_array("2nd", $positions)) { $p_pos1 = "1";$p_pos3 = "0";}
					elseif ( ! in_array("3rd", $positions)) { $p_pos1 = "1";$p_pos2 = "0";}
				}
			}
			else { 
				if (in_array("all", $positions)) {
					$p_pos1 = "0";$p_pos2 = "1";$p_pos3 = "2";
				}
				else { 
					if ( ! in_array("1st", $positions)) { $p_pos2 = "0";$p_pos3 = "1";}
					elseif ( ! in_array("2nd", $positions)) { $p_pos1 = "0";$p_pos3 = "1";}
					elseif ( ! in_array("3rd", $positions)) { $p_pos1 = "0";$p_pos2 = "1";}
				}
			}
			// set end of gene
			$phybp_end = $phybp + $charset_count[$gCPHY] - 1; 
			// always includes "asone" GENE partitions for NEXUS (with name of specified codon position/s)
			if ( in_array("all", $positions) || in_array("aas", $positions)){
				$nex_partitions[] = "\tcharset $gCPHY = $phybp-$phybp_end;";
				$nex_all_partitions[] = $gCPHY;
				if ($by_positions != "each" && $by_positions !="123") {
					$which_nex_partitions[] = $gCPHY;
				}
			}
			else{
				$nex_partitions[] = "\tcharset $gCPHY" . "_pos$part_name_end" . " = $phybp-$phybp_end;";
				$nex_all_partitions[] = $gCPHY . "_pos$part_name_end";
				if ($by_positions != "each" && $by_positions !="123") {
					$which_nex_partitions[] = $gCPHY . "_pos$part_name_end";
				}
			}
			// and include "asone" GENE partition for special
			// if ( isset($gene_by_positions) && isset($gene_positions) && $by_positions == "asone"){
				// if ( in_array("all", $positions)) {$which_nex_partitions[] = $gCPHY;}
				// else {
					//$nex_partitions[] = "\tcharset $gCPHY" . "_pos$part_name_end" . " = $phybp-$phybp_end;";
					// $which_nex_partitions[] = $gCPHY . "_pos$part_name_end";
				// }
			// }
			if ($by_positions == "asone" ){ //for simple gene partitions
				if (in_array("all", $positions)){
					$phy_partitions[] = "DNA, $gCPHY = $phybp-$phybp_end";
				}
				elseif (in_array("aas", $positions)){
					if($aa_or_not[$gCPHY] == 'yes'){
						$phy_partitions[] = "JJT, $gCPHY = $phybp-$phybp_end";
					}
					else {$phy_partitions[] = "DNA, $gCPHY = $phybp-$phybp_end";}
				}
				else{
					$phy_partitions[] = "DNA, $gCPHY". "_pos$part_name_end = $phybp-$phybp_end";
				}
			}
			elseif ( $by_positions == "each" ){ // for single codon positon partitions
				foreach ( $positions  as $act_pos ) {
					if ( $act_pos == "1st" || $act_pos == "all" ){
						$pos1_start = $phybp + $p_pos1;
						$phy_partitions[] = "DNA, " . $gCPHY . "_pos1 = " . $pos1_start . "-" . $phybp_end .  $p_jumps ;
						$nex_partitions[] = "\tcharset " . $gCPHY . "_pos1 = " . $pos1_start . "-" . $phybp_end .  $p_jumps . ";" ;
						$which_nex_partitions[] = $gCPHY . "_pos1";
					}
					if ( $act_pos == "all" || $act_pos == "2nd" ){
						$pos2_start = $phybp + $p_pos2;
						$phy_partitions[] = "DNA, " . $gCPHY . "_pos2 = " . $pos2_start . "-" . $phybp_end .  $p_jumps ;
						$nex_partitions[] = "\tcharset " . $gCPHY . "_pos2 = " . $pos2_start . "-" . $phybp_end .  $p_jumps . ";" ;
						$which_nex_partitions[] = $gCPHY . "_pos2";
					}
					if ( $act_pos == "all" || $act_pos == "3rd" ){
						$pos3_start = $phybp + $p_pos3;
						$phy_partitions[] = "DNA, " . $gCPHY . "_pos3 = " . $pos3_start . "-" . $phybp_end .  $p_jumps ;
						$nex_partitions[] = "\tcharset " . $gCPHY . "_pos3 = " . $pos3_start . "-" . $phybp_end .  $p_jumps . ";" ;
						$which_nex_partitions[] = $gCPHY . "_pos3";
					}
					}
				}
			else { //partion 1 with codon 1+2 and partition2 with codon 3
				$pos1_start = $phybp + $p_pos1; $pos2_start = $phybp + $p_pos2; $pos3_start = $phybp + $p_pos3;
				$phy_partitions[] = "DNA, " . $gCPHY . "_pos12 = " . $pos1_start . "-" . $phybp_end . "\\3, " . $pos2_start . "-" . $phybp_end . "\\3"; 
				$phy_partitions[] = "DNA, " . $gCPHY . "_pos3 = " . $pos3_start . "-" . $phybp_end . "\\3";
				$which_nex_partitions[] = $gCPHY . "_pos12";
					if ( $pos1_start < $pos2_start ) { $pos12_start = $pos2_start; } else { $pos12_start = $pos1_start; }
				$nex_partitions[] = "\tcharset " . $gCPHY . "_pos12 = " . $pos1_start . "-" . $phybp_end . "\\3 " . $pos2_start . "-" . $phybp_end . "\\3;";
				$nex_partitions[] = "\tcharset " . $gCPHY . "_pos3 = " . $pos3_start . "-" . $phybp_end . "\\3;";
				$which_nex_partitions[] = $gCPHY . "_pos3";
			}
			$phybp = $phybp + $charset_count[$gCPHY];
			//adding intron partitions
			if (count($intron_dataset[$gCPHY]) > 0 && $ignore_introns == 'no'){
				$intr = 1;
				foreach ($intron_lengths[$gCPHY] as $il) {
					$nex_partitions[] = "\tcharset " . $gCPHY . "_intron_$intr = " . $phybp . "-" . ($phybp+$il-1) . ";" ;
					$which_nex_partitions[] = $gCPHY . "_intron_$intr";
					$nex_all_partitions[] = $gCPHY . "_intron_$intr";
					$phy_partitions[] = "DNA, " . $gCPHY . "_intron_$intr = " . $phybp . "-" . ($phybp+$il-1);
					$phybp = $phybp+$il;
					$intr++;
					
				}
			}
		}
	}

	// #################################################################################
	// Section: Build output - end lines and end info
	// #################################################################################
	// creating list of partitions to enter into appropriate places in nex output
	if( $format == "NEXUS" ) {
		// set counts on partitions
		$count_nex = 0;
		$part_list_sum = array();
		$i = 0;
		foreach ($geneCodes as $gCPHY){
			$i = $i+1;
			//echo "$gCPHY = ".$aa_or_not[$gCPHY]."</br>";
			if ($aa_or_not[$gCPHY] == 'yes'){ $part_list_sum['full_aa'][] = $i;}
			else { $part_list_sum['full_dna'][] = $i;}
			if (count($intron_dataset[$gCPHY]) > 0 && $ignore_introns == 'no'){
				foreach ($intron_lengths[$gCPHY] as $ils) {
					$i= $i+1;
					$part_list_sum['full_dna'][] = $i;
				}
			}
		}
		$i = 0;
		foreach ($which_nex_partitions as $wnp){
			$i = $i+1;
			if ($pos = strpos($wnp, "_pos") == FALSE) {
				if ($aa_or_not[$wnp] == 'yes'){
					$part_list_sum['special_aa'][] = $i;
				}
				else { $part_list_sum['special_dna'][] = $i;}
			}
			else { $part_list_sum['special_dna'][] = $i;}
		}
		//print_r ($part_list_sum);echo "</br>";
		if (count($part_list_sum[full_dna])==0){$part_list_sum[full_dna][] = "";} 
	}
	// creating end block of file
	if( $format == "TNT" ) { 
		$output .= ";\nproc/;";
	}
	elseif( $format == "FASTA" ) {
		$output .= "";
	}
	elseif( $format == "PHYLIP" ) {
		$output .= "";
		$phy_partitions = implode("\n", $phy_partitions);
	}
	else {
		
		$output .= ";\nEND;\n\n";
		
		$output .= "begin mrbayes;\n";
		
		$i = 1;
		$a = 0;
		$b = 0;
		$output .= implode("\n", $nex_partitions);
		$output_genes = implode(", ", $geneCodes);
		$output_genes = implode(", ", $nex_all_partitions);
		$output .= "\npartition GENES = " . count($nex_all_partitions). ": $output_genes;";
		
		if ( $by_positions != "asone" && !isset($gene_by_positions) && !isset($gene_positions)){
			$output_parts = implode(", ", $which_nex_partitions);
			$output .= "\npartition CODONPOSITIONS = " . count($which_nex_partitions) . ": $output_parts;";
			$output .= "\n\nset partition = CODONPOSITIONS;\n";
		}
		elseif ( isset($gene_by_positions) && isset($gene_positions) && $which_nex_partitions != $nex_all_partitions){
			$output_parts = implode(", ", $which_nex_partitions);
			$output .= "\npartition SPECIAL = " . count($which_nex_partitions) . ": $output_parts;";
			$output .= "\n\nset partition = SPECIAL;\n";
		}
		else {
			$output .= "\n\nset partition = GENES;\n";
		}
		$output .= "\nset autoclose=yes;\n";
		if( isset($outgroup) ) { 
			$output .= "outgroup $outgroup;\n"; 
		}
		$output .= "prset applyto=(all) ratepr=variable brlensp=unconstrained:Exp(100.0) ";
		$output .= "shapepr=exp(1.0)";
		if ($part_list_sum['full_dna'][0] != "" && $part_list_sum['special_dna'] > 0){
			$output .= " tratiopr=beta(2.0,1.0)";
		}
		$output .= ";\n";
		if( in_array("yes", $aa_or_not) || in_array("aas", $positions) ) {
		//print_r ($aa_or_not);
			if( $which_nex_partitions != $nex_all_partitions ) {
				$output .= "prset applyto=(". implode(',', $part_list_sum['special_aa'])  ." [". implode(',', $part_list_sum['full_aa'])  ."]) aamodelpr=mixed;\n";
				$output .= "lset  applyto=(". implode(',', $part_list_sum['special_aa'])  ." [". implode(',', $part_list_sum['full_aa'])  ."]) rates=gamma [invgamma];\n";
				if (in_array("no", $aa_or_not) || $ignore_introns == 'no' && count($intron_lengths) > 0){
					$output .= "lset  applyto=(". implode(',', $part_list_sum['special_dna']) ." [". implode(',', $part_list_sum['full_dna']) ."]) nst=mixed rates=gamma [invgamma];\n";
					$output .= "unlink statefreq=(". implode(',', $part_list_sum['special_dna']) ." [". implode(',', $part_list_sum['full_dna']) ."]);\n";
					}
			}
			else { 
				$output .= "prset applyto=(". implode(',', $part_list_sum['full_aa']) .") aamodelpr=mixed;\n";
				$output .= "lset  applyto=(". implode(',', $part_list_sum['full_aa'])  .") rates=gamma [invgamma];\n";
				if (in_array("no", $aa_or_not) || $ignore_introns == 'no' && count($intron_lengths) > 0 ){
					$output .= "lset  applyto=(". implode(',', $part_list_sum['full_dna']) .") nst=mixed rates=gamma [invgamma];\n";
					$output .= "unlink statefreq=(". implode(',', $part_list_sum['full_dna']) .");\n";
				}
			}
		}
		else {
			$output .= "lset applyto=(all) nst=mixed rates=gamma [invgamma];\n";
			$output .= "unlink statefreq=(all);\n";
		}
		
		$output .= "unlink shape=(all) revmat=(all) tratio=(all) [pinvar=(all)];\n";
		$output .= "mcmc ngen=10000000 printfreq=1000 samplefreq=1000 nchains=4 nruns=2 savebrlens=yes [temp=0.11];\n";	
        $output .= " sump relburnin=yes [no]  burninfrac=0.25 [2500];\n";
        $output .= " sumt relburnin=yes [no]  burninfrac=0.25 [2500] contype=halfcompat [allcompat];\n";
        $output .= "end;";
	}

	// #################################################################################
	// Section: HTML output and download links
	// #################################################################################
	// print html headers
	$admin = false;
	$in_includes = true;
	include_once 'header.php';
	// print navigation bar
	nav();
	$output2 = $output;
	if ( $format == "PHYLIP" && count($geneCodes) < 2 && $by_positions == "asone") { unset($phy_partitions); }


	// begin HTML page content
	echo "<div id=\"content\">";

	if( count($warning) > 0 ) {
		// issue warnings using javascript
		echo "<script type=\"text/javascript\">\n";
		foreach( $warning as $item ) {
			echo "alert(\"" . $item . "\")\n";
		}
		echo "</script>";
	}
	?>
			<form action="dataset_to_file.php" method="post">
			<table border="0" width="960px" cellpadding="5px"> <!-- super table -->
			<tr>
				<td align='center' width="100%" class="label4">Your dataset: copy/edit/create file</td>
			</tr>
			<tr>
				<td class="field1"><textarea rows="35" cols="125" wrap='off' name="dataset"><?php echo $output2; ?></textarea></td>
			</tr>
				<tr><input type="hidden" name="format" value="<?php echo $format; ?>">
			<td>
			<input type="submit" name="submit" value="Create file">
			<?php if ( $format == "PHYLIP" && isset($phy_partitions)) {
				?><!--<input type="hidden" name="phy_partitions" value="< ?php echo $phy_partitions; ?>"> -->
				<tr>
				<td class="field"><textarea rows="10" cols="40" wrap='off' name="phy_partitions" ><?php echo $phy_partitions; ?></textarea></td>
				</tr><tr>
				<td><input type="submit" name="phy_parts" value="Create PHYLIP partitions file"><?php } ?>
			</td>
		</tr>
	</table>
			</form>
	</div> <!-- end content -->

	<?php
	make_footer($date_timezone, $config_sitename, $version, $base_url);

	?>
	</body>
	</html>
	<?php
	}
	}
}
?>
