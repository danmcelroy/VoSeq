<?php
// #################################################################################
// #################################################################################
// Voseq includes/make_fasta_genbank.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Creates a fasta list with voucher info and sequence
// for GenBank submission
// #################################################################################
// #################################################################################
// Section: Startup/includes
// #################################################################################
// includes
ob_start();//Hook output buffer - disallows web printing of file info...
include'../conf.php';
ob_end_clean();//Clear output buffer//includes
include '../functions.php';
include 'translation_functions.php';
//include'../markup-functions.php';
// #################################################################################
// Section: Functions - clean_item() and show_errors()
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
function show_errors($se_in) {
		// error found
			include'../markup-functions.php';
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

function do_alert($msg) 
{
	echo '<script type="text/javascript">alert("' . $msg . '"); </script>';
}


// #################################################################################
// Section: Function - process_fasta_sequence
// #################################################################################
// @brief: takes a sequence string and puts it in frame (starts with 1st codon position
// @input:  string sequence
// @output: string sequence replaces ?s with Ns
//			remove ?s from the end
// TODO: put it in frame
function process_fasta_sequence($sequences) {
	$sequences = str_replace("?", "N", $sequences);
	$sequences = preg_replace("/N+$/", "", $sequences);
	return $sequences;
}



// open database connections
// #################################################################################
// Section: Get code(s) and gene(s) and taxonomic info
// #################################################################################
// open database connections
@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect');
mysql_select_db($db) or die ('Unable to select database');
if( function_exists(mysql_set_charset) ) {
	mysql_set_charset("utf8");
}


if (trim($_POST['codes']) != ""){
	$raw_codes = explode("\n", $_POST['codes']);
}else{ unset($raw_codes); }
//$raw_codes = split("\n", $_POST['codes']);

// geneCodes here
unset($geneCodes);
if (isset($_POST['geneCodes'])){
	foreach ( $_POST['geneCodes'] as $k1=> $c1){ //putting choosen genes into array
		if ($c1 == 'on')	{
			$genes[] =  $k1;
		}
	}
}elseif ($_POST['genesets'] == "Choose geneset"){
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
if (isset($geneset_taxa) && isset($genes)){$genes = array_merge( $geneset_taxa, $genes) ;}
elseif (isset($geneset_taxa) && ! isset($genes)){$genes = $geneset_taxa ;}
elseif (! isset($geneset_taxa) && isset($genes)){$genes = $genes ;}
else { $errorList[] = "No genes are chosen!</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Pointless to make a table without genes..."; }


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

// codes here
$lines = array();

if (isset($raw_codes)){
$raw_codes = array_unique($raw_codes);
foreach($raw_codes AS $item) {
	$item = clean_item($item);
	$item = trim($item);
	if ($item != "") {
		$cquery = "SELECT code FROM ". $p_ . "vouchers WHERE code='$item'";
		$cresult = mysql_query($cquery) or die("Error in query: $query. " . mysql_error());
		// if records present
		if( mysql_num_rows($cresult) > 0 ) {
			while( $row = mysql_fetch_object($cresult) ) {		
				array_push($lines, $item);
			}
		}
		else {
		$errorList[] = "No voucher named <b>$item</b> exists in database!</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please add it in the voucher section or remove it from taxon set!";
		}
	}
}unset($item);

$lines = array_unique($lines);
}
//check for error and if none proceed with building table
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
// Section: Build fasta list
// #################################################################################
$output = ""; // for sequence file
$protout = ""; // for protein file if proteincoding and introns
$numacc = array(); //counting number of taxa with accession number
foreach($genes as $geneCode) {
	if ($output != "") {
		$output .= "\n\n";
	}
	if ($protout != "") {
		$protout .= "\n\n";
	}
	foreach($lines as $code) { // iterate through codes
		$query = "SELECT " . $p_ . "vouchers.orden, 
						 " . $p_ . "vouchers.code,
						 " . $p_ . "vouchers.family, 
						 " . $p_ . "vouchers.subfamily, 
						 " . $p_ . "vouchers.tribe, 
						 " . $p_ . "vouchers.genus, 
						 " . $p_ . "vouchers.species, 
						 " . $p_ . "sequences.sequences, 
						 " . $p_ . "sequences.accession, 
						 " . $p_ . "genes.readingframe, 
						 " . $p_ . "genes.length,
						 " . $p_ . "genes.genetic_code,
						 " . $p_ . "genes.intron, 
						 " . $p_ . "genes.description FROM 
						 " . $p_ . "vouchers, 
						 " . $p_ . "sequences, 
						 " . $p_ . "genes WHERE 
						 " . $p_ . "vouchers.code='$code' AND 
						 " . $p_ . "sequences.code='$code' AND 
						 " . $p_ . "sequences.geneCode='$geneCode' AND 
						 " . $p_ . "genes.geneCode='$geneCode'";
		$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());

		if (mysql_num_rows($result) > 0) { // get results from MySQL
			unset($lineage); 
			$lineage = " [Lineage=";
			while ($row = mysql_fetch_object($result)) {
				if ($row->accession !== "" && trim($row->accession) !== "" && $row->accession !== "NULL"){
					$numacc[] = $code . " for " . $geneCode . " = " . $row->accession ;}
				else {
					if( $row->orden ) {
						$lineage .= " $row->orden;";
					}
					if( $row->family ) {
						$lineage .= " $row->family;";
					}
					if( $row->subfamily ) {
						$lineage .= " $row->subfamily;";
					}
					if( $row->tribe ) {
						$lineage .= " $row->tribe;";
					}
					$lineage .= " $row->genus] ";
					$species = str_replace(" ", "_", $row->species);
					$outputS = ">" . $row->genus . "_" . $species . "_";
					$outputS .= $row->code . " [org=$row->genus $row->species] ";
					$outputS .= "[Specimen-voucher=$row->code]";
					$outputS .= " [note=" . $row->description . " gene, partial cds.]";
					$outputS .= " $lineage";
					// need to replace ? with N and put it in frame 
					// (sequence starts with 1st codon position)
					$sequences = $row->sequences;
					$sequences = process_fasta_sequence($sequences);
					if ($row->intron == "" || $row->intron == "NULL"){
						$outputS .= "\n$sequences\n";
					}
					else {
						$outputP = ">" . $row->genus . "_" . $species . "_";
						$outputP .= $row->code . " [gene=$geneCode] [protein=$row->description] ";
						#intron array
						if (strlen($sequences) < $row->length){
							$sequences = str_pad($sequences, $row->length, "?");
						}
						// translating for prootein file
						$prot = remove_introns($row->sequences,$row->intron);
						if ($row->readingframe == "2") { $prot = substr($prot[0],1);}
						elseif ($row->readingframe == "3") { $prot = substr($prot[0],2);}
						else {$prot = $prot[0];}
						$tprot=translate_DNA_to_protein($prot,$row->genetic_code);//translate
						$clean_tprot = preg_replace(array("/^X+/","/X+$/"),"",$tprot);
						$outputP .= "\n" . $clean_tprot ."\n";
						// removing introns and replacing chars for seq file
						$introns = explode(";", $row->intron);
						$currpos = 0;
						$new_seq_arr = array();
						foreach ($introns as $intr){
																		#$new_seq_arr[] = "\n".$intr."\n";
							$ise = explode("-", $intr);
							$ise[2] = $ise[1] - $ise[0] + 1;
							if ($ise[0] !== "1"){
																		#$new_seq_arr[] = "Exon:";
																		#$new_seq_arr[] = "\n" . $currpos ." v " . ($ise[0]-1-$currpos) ."\n";
								$new_seq_arr[] = str_replace(array("-","?"),"N",substr($sequences,$currpos,$ise[0]-1-$currpos));
							}
																		#$new_seq_arr[] = "\nIntron:\n" . ($ise[0]-1) ." w ".$ise[2] ."\n";
							$new_seq_arr[] = str_replace(array("N","?"),"-",substr($sequences,$ise[0]-1,$ise[2]));
							$currpos = $ise[1];
							if ($intr == end($introns)){
																		#$new_seq_arr[] = "\nsista Exon\n" . $currpos ." x " . $row->length ."\n";
								$new_seq_arr[] = str_replace(array("-","?"),"N",substr($sequences,$currpos,$row->length-$currpos));
							}
						}
						$sequences = implode("",$new_seq_arr);
						$outputS .= "\n$sequences\n";
					}
					$output .= $outputS;
					$protout .= $outputP;
				}
			}
		}
	}
}

// #################################################################################
// Section: Show output or error message
// #################################################################################
if ($output == "" ){
	$errorList[] = "No voucher had any of the choosen genes!</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please choose existing voucher-gene combinations...";
	$title = "$config_sitename: Dataset Error";
	// print html headers
	$admin = false;
	$in_includes = true;
	include_once 'header.php';
	//print errors
	show_errors($errorList);
}
else{ //start building dataset
	if( $php_version == "5" ) {
		// date_default_timezone_set($date_timezone); php5
		date_default_timezone_set($date_timezone);
	}

	# filename for download prot file
	//date_default_timezone_set($date_timezone);php5
	// if ($protout != ""){
		// $genbank_file = "genbank_file_" . date('Ymd') . ".txt"; 
		// header("Content-Type: application/vnd.ms-notepad");
		// header("Content-Disposition: attachment; filename=$genbank_file");
		// echo "[You have choosen to make GenBank submission in Sequin with intron data.\nBelow is first a part including sequence fasta lines and further 
			// on a part with protein data.\n Cut out the protein part and put in another text file for use in Sequin.]\n\n\n\n$output\n\n\n\n$protout";
	// }
	// else{
		// $genbank_file = "genbank_file_" . date('Ymd') . ".txt"; 
		// header("Content-Type: application/vnd.ms-notepad");
		// header("Content-Disposition: attachment; filename=$genbank_file");
		// echo $output;
	// }
		// //warning bout number of already accessionized taxa
		// ? >< script type="text/javascript">
			// ("taxa already have accession numbers. \nThese will not be prepared here!");
		// <script><?php
	// #################################################################################
	// Section: HTML output and download links
	// #################################################################################
	// print html headers
	$admin = false;
	$in_includes = true;
	include_once 'header.php';
	include "../markup-functions.php";
	// print navigation bar
	nav();
	$output2 = str_replace("\x1F","",$output);
	if ( $format == "PHYLIP" && count($geneCodes) < 2 && $by_positions == "asone") { unset($phy_partitions); }

	// begin HTML page content
	echo "<div id=\"content\">";
	?>
			 <!-- super table -->
			<?php if( count($numacc) > 0 ) {
				// issue warnings using javascript
				$item = implode("\n",$numacc);
				$noitem = count($numacc);
				//echo "$noitem\n$item";
				echo "<table width='100px' cellpadding='5px' align='center'><tr><td align='center' width='10%' class='label4'>$noitem sequences already have accession.<br>These will be ignored:</td></tr><tr>
				<td class='field1' align='center'><textarea rows='10' cols='50' wrap='off' name='acced'>" . $item . "</textarea><br></td></tr></table>";
			} ?>
			<form action="dataset_to_file.php" method="post">
			<table border="0" width="960px" cellpadding="5px">
			<td><h1>This is your nucleotide fasta file to import into Sequin (http://www.ncbi.nlm.nih.gov/Sequin/).</h1>
			<input type="submit" name="submit_gB" value="Create nucleotide file"></td>
			<tr>
				<td align='center' width="100%" class="label4">Your nucleotide file:</td>
			</tr>
			<tr>
				<td class="field1"><textarea rows="30" cols="125" wrap='off' name="output"><?php echo $output; ?></textarea></td>
			</tr>
				<tr><input type="hidden" name="format" value="<?php echo $format; ?>">
			
			<?php if ($protout != ""){
				?><tr>
				<td><h1>This is your protein file to use in case you have introns. </br>
				Load it in Sequin after your nucleotide file to link nucleotides to protein and identify the introns.</h1>
				<input type="submit" name="prot" value="Create protein file"><td>
				</tr>
				<tr>
					<td align='center' width="100%" class="label4">Your protein file:</td>
				</tr>
				<tr>
					<td class="field"><textarea rows="20" cols="125" wrap='off' name="protout" ><?php echo $protout; ?></textarea></td>
				</tr><tr>
				<?php } ?>
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
?>
