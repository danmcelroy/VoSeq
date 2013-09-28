<?php
// ############################################################################
// ############################################################################
// Voseq admin/mysqlmerger.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: process merger of two databases and create result output.
// Connected to mysqlmerge.php
//
// ############################################################################


// #################################################################################
// Section: Create downloadable file with table
// #################################################################################

// Name of the file
$filename = "uploads/" . trim($_POST['name']);
$filename = preg_replace('/[^\w\._^\/]+/', '_', $filename);

if( !file_exists($filename) ) {
	echo "Error, couldn't find that file!";
	die(0);
}

//check admin login session
include'../login/auth-admin.php';
// includes
ob_start();//Hook output buffer - disallows web printing of file info...
include'../conf.php';
ob_end_clean();//Clear output buffer//includes

// set increased time limit for larger calculations -> 1000 seconds
set_time_limit(1000);


// Connect to MySQL server
mysql_connect($host, $user, $pass) or die('Error connecting to MySQL server: ' . mysql_error());
// Select database
mysql_select_db($db) or die('Error selecting MySQL database: ' . mysql_error());

$tablearraynames = array("members", "vouchers", "genes", "sequences", "primers", "taxonsets", "genesets","search","search_results");
foreach ($tablearraynames as $taname){
	
	$q_droptemps = "DROP TABLE IF EXISTS temp_".$p_.$taname.";";
	$r_droptemps = mysql_query($q_droptemps) or die ("Error in query(11): $q_droptemps. " . mysql_error());
}
// Temporary variable, used to store current query
$templine = '';
// Read in entire file
$lines = file($filename);
// array for collection of table names
$tablenames = array();
$tests=array();
// Loop through each line
foreach ($lines as $line)
{
	// Skip it if it's a comment
	if (substr($line, 0, 2) == '--' || $line == '' || substr($line, 0, 3) == '/*!') {
		continue;
	}
	// changing table names to temporary import tables
	if (strpos($line,'DROP TABLE IF EXISTS') !== FALSE || strpos($line,'CREATE TABLE') !== FALSE || strpos($line,'LOCK TABLES') !== FALSE || strpos($line, 'INSERT INTO') !== FALSE){
		preg_match('/\`[^\`]*\`/',$line, $tn);
		//echo "<br><br>Drop table match:" . $tn[0]."<br>";
		$tn = str_replace('`','',$tn[0]);
		$tablenames[] = $tn;
		if (strpos($tn, "members") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "members", $line); }
		if (strpos($tn, "vouchers") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "vouchers", $line); }
		if (strpos($tn, "genesets") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "genesets", $line); }
		elseif (strpos($tn, "genes") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "genes", $line); }
		if (strpos($tn, "sequences") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "sequences", $line); }
		if (strpos($tn, "primers") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "primers", $line); }
		if (strpos($tn, "taxonsets") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "taxonsets", $line); }
		if (strpos($tn, "search_results") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "search_results", $line); }
		elseif (strpos($tn, "search") !== FALSE) {$line = str_replace($tn, "temp_" . $p_ . "search", $line); }
		//echo substr($line,0,50) . "\n";
	}
	// Add this line to the current segment
	$templine .= $line;
	// If it has a semicolon at the end, it's the end of the query
	if (substr(trim($line), -1, 1) == ';')
	{
		// Perform the query
		 mysql_query($templine) or die ('Error performing query(1) \'<strong>' . $templine . '\': ' . mysql_error() . '<br /><br />');
		 //or print('Error: tnames:'. print_r($tablenames) . "<br>tests: ".print_r($tests).'<br><br>');
												//
		// Reset temp variable to empty
		$templine = '';
	}
}

//function find new uniques
function find_news ($table, $main_headers, $comparelist) {
	$mainheaders = implode(",",$main_headers);
	$q1 = "SELECT $mainheaders FROM temp_$table";
	//echo"<br>$q1<br>";
	$r1 = mysql_query($q1) or die("Error in query(2a): $q1. " . mysql_error());
	// if records present
	$new_mainheaders = array();
	if( mysql_num_rows($r1) > 0 ) {
		while( $row1 = mysql_fetch_object($r1) ) {
			if (count($main_headers) > 1){
				$new_mainheaders[] = array($row1->$main_headers[0],$row1->$main_headers[1]) ;
			}
			else {$new_mainheaders[] = $row1->$main_headers[0];}
		}
	}
	//print_r($new_mainheaders);
	// find new mainheaders (= not present in old ones), including second mainheader
	$new_mainheaders2 = array();
	$mod_mainheaders2 = array();
	$update_mainheaders2 = array();
	$same_mainheaders2 = array();
	foreach ($new_mainheaders as $nmh1) {
		if (count($main_headers) > 1){ 
			$q2 = "SELECT * FROM $table WHERE ". $main_headers[0] . "='".$nmh1[0]."' AND ". $main_headers[1] ."='".$nmh1[1]."'";
			$q2temp = "SELECT * FROM temp_$table WHERE ". $main_headers[0] . "='".$nmh1[0]."' AND ". $main_headers[1] ."='".$nmh1[1]."'";
		}
		else{ 
			$q2 = "SELECT * FROM $table WHERE ". $main_headers[0] . "='$nmh1'";
			$q2temp = "SELECT * FROM temp_$table WHERE ". $main_headers[0] . "='$nmh1'";
		}
		//echo "<br>$q2<br>$q2temp<br>";
		$r2 = mysql_query($q2) or die("Error in query(2b): $q2. " . mysql_error());
		$r2temp = mysql_query($q2temp) or die("Error in query(3): $q2temp. " . mysql_error());
		if( mysql_num_rows($r2) > 0 ) {
		//if records present
			while( $row2 = mysql_fetch_object($r2) ) {
				while( $row2temp = mysql_fetch_object($r2temp) ) {
				$same = 0;
					foreach ($comparelist as $cl){
						if ($row2temp->$cl != "" && $row2->$cl != ""&& $row2->$cl != $row2temp->$cl){
							//echo "<br>".$nmh1." ($cl) = ".$row2->$cl." + ". $row2temp->$cl."<br>";
							$same = 1;
							if (count($main_headers) > 1){
								$mod_mainheaders2[$nmh1[0]."§".$nmh1[1]][] = $cl;
							}
							else{
								$mod_mainheaders2[$nmh1][] = $cl;
							}
						}
						elseif ($row2->$cl == "" && $row2temp->$cl != "" && $row2->$cl != $row2temp->$cl){
							$same = 1;
							//echo "<br>".$nmh1." ($cl) = ".$row2->$cl." + ". $row2temp->$cl."<br>";
							if (count($main_headers) > 1){
								$update_mainheaders2[$nmh1[0]."§".$nmh1[1]][] = $cl;
							}
							else{
								$update_mainheaders2[$nmh1][] = $cl;
							}
						}
						else{//echo "<br>".$nmh1." ($cl) = ".$row2->$cl." + ". $row2temp->$cl."<br>";
						}
					}
					if ($same == 0){
						if (count($main_headers) > 1){
							$same_mainheaders2[] = $nmh1[0]."§".$nmh1[1];
						}
						else{
							$same_mainheaders2[$nmh1] = "X";
						}
					}
				}
			}
		}
		else{
			if (count($main_headers) > 1){
				$new_mainheaders2[] = $nmh1[0]."§".$nmh1[1];
			}
			else{
				$new_mainheaders2[] = $nmh1;
			}
		}
	}
	$same_mainheaders2 = count($same_mainheaders2);
	return array($new_mainheaders2,$update_mainheaders2,$mod_mainheaders2,$same_mainheaders2);//, $moddednames);
}

// get new values using function: find_news

 $members_news = find_news($p_.'members', array('login'), array('firstname','lastname','passwd', 'admin'));
 $genes_news = find_news($p_.'genes', array('geneCode'), array('length','readingframe','genetype','prot_code','intron','aligned','genetic_code'));
 $sequences_news = find_news($p_.'sequences', array('code','geneCode'), array('labPerson','sequences','notes','accession','genbank'));
 $primers_news = find_news($p_.'primers', array('code','geneCode'), array('primer1','primer2','primer3','primer4','primer5','primer6'));
 $vouchers_news = find_news($p_.'vouchers', array('code'), array('orden','family','subfamily','genus','species','subspecies','country','specificLocality','latitude','longitude','altitude','collector','dateCollection','determinedBy','voucherImage','thumbnail','extraction','dateExtraction','extractor','voucherLocality','publishedIn','notes','hostorg','sex','extractionTube','voucher','voucherCode','flickr_id','auctor'));
 $taxonsets_news = find_news($p_.'taxonsets', array('taxonset_name'), array( 'taxonset_list','taxonset_creator','taxonset_description'));
 $genesets_news = find_news($p_.'genesets', array('geneset_name'), array( 'geneset_list','geneset_creator','geneset_description'));


// setting up bad value arrays (of vouchers and genes) that will affect others. Bad values are vouchers and genes that are not added due to conflicts
$bad_vouchers = array();
if (count($vouchers_news[2])>0){
	foreach ($vouchers_news[2] as $vn2 => $vn2h){
			$bad_vouchers[] = $vn2;
	}
}
$bad_genes = array();
if (count($genes_news[2])>0){
	foreach ($genes_news[2] as $vn2 => $vn2h){
			$bad_genes[] = $vn2;
	}
}
// start outputs and handling variables
$output = '<b><u>This is an output from the merging of database1 with database2</b></u><br>';
$output_bads = '<br><b>These sequences/primers are not uploaded due to problems with connected codes/genes:</b><br>';
$output0 = '<br><b>These records are added as new:</b><br>';
$output1 = '<br><b>These records are updated with added information:</b><br>';
$output2 = '<br><b>These records are not updated due to conflict with existing records:</b><br>';
$output3 = '<br><b>These records are ignored as they already exist:</b><br>';
$tablearraynames = array("members", "vouchers", "genes", "sequences", "primers", "taxonsets", "genesets");
$tablearrayvitals = array("login", "code", "geneCode", array("code","geneCode"), array("code","geneCode"), "taxonset_name", "geneset_name");
$tablearray = array($members_news, $vouchers_news, $genes_news, $sequences_news, $primers_news, $taxonsets_news, $genesets_news);


for ($i = 0; $i <= 6; $i++){
	$newsarr = $tablearray[$i];
	//echo"<br><br>ACTIVE ARR IS:"; echo strtoupper($tablearraynames[$i])."<br><br>";
	//checking bad data and moves sequence or primer records from $newsarr[0] and $newarr[1] to a new $newsarr[4]
	if ($tablearraynames[$i] == 'sequences' || $tablearraynames[$i] == 'primers' ){
		//get existing vouchers
		$btv = array();
		$qbadtestv = "SELECT code FROM ".$p_."vouchers";
		$rbadtestv = mysql_query($qbadtestv) or die("Error in query(4): $qbadtestv. " . mysql_error());
		//if records present
		if( mysql_num_rows($rbadtestv) > 0 ) {
			while( $rowbadtestv = mysql_fetch_object($rbadtestv) ) {
				foreach ($rowbadtestv as $rbtv){
					$btv[] = $rbtv;
				}
			}
		}
		//get existing genes
		$btg = array();
		$qbadtestg = "SELECT geneCode FROM ".$p_."genes";
		$rbadtestg = mysql_query($qbadtestg) or die("Error in query(5): $qbadtestg. " . mysql_error());
		//if records present
		if( mysql_num_rows($rbadtestg) > 0 ) {
			while( $rowbadtestg = mysql_fetch_object($rbadtestg) ) {
				foreach ($rowbadtestg as $rbtg){
					$btg[] = $rbtg;
				}
			}
		}
		foreach ($newsarr[0] as $nwannum => $nwan){
			$tav = explode("§",$nwan);
			//echo $tav[0]." ".$tav[1]."<br>";
			if (in_array($tav[0], $bad_vouchers) || in_array($tav[1], $bad_genes) || !in_array($tav[0],$btv) || !in_array($tav[1],$btg)){
				$output_bads .= $tablearraynames[$i] ." for ". $tav[0] ." and ". $tav[1] ."<br>";
				unset($newsarr[0][$nwannum]);
			}
		}
		//echo "<br>1<br>";
		//print_r($btg);print_r($btv);
		foreach ($newsarr[1] as $nwan => $nwanh){
			$tav = explode("§",$nwan);
			//echo $tav[0]." ".$tav[1]."<br>";
			if (in_array($tav[0], $bad_vouchers) || in_array($tav[1], $bad_genes) || !in_array($tav[0],$btv) || !in_array($tav[1],$btg)){
				$output_bads .= $tablearraynames[$i] ." for ". $tav[0] ." and ". $tav[1] ."<br>";
				unset($newsarr[1][$nwan]);
			}
		}
	}
	// if (count($vouchers_news[0])>0 || count($vouchers_news[1])>0 ||count($vouchers_news[2])>0 ||count($vouchers_news[3])>0){
		// $output .= "<br><br>From " . strtoupper($tablearraynames[$i]) . ":<br>";
	// }
	//Start inputting and outputting
	// going through and inputs new records
	if (count($newsarr[0])>0){
		$output0 .= "<br>New added ".$tablearraynames[$i].":<br>";
		//getting columns that are shared between databases
		$cols = array();
		$vtable_array = array($p_.$tablearraynames[$i],"temp_".$p_.$tablearraynames[$i]);
		foreach ($vtable_array as $vv){
			$qc = "SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'temp_".$p_.$tablearraynames[$i]."'";
			$rc = mysql_query($qc) or die("Error in query(6): $qc. " . mysql_error());
			//if records present
			if( mysql_num_rows($rc) > 0 ) {
				while( $rowc = mysql_fetch_object($rc) ) {
					foreach ($rowc as $rch){
						$cols[$vv][$rch] = $rch;
					}
				}
			}
		}
		$aint = array_intersect($cols[$vtable_array[1]],$cols[$vtable_array[0]]);
		unset ($aint['id'],$aint['member_id'],$aint['taxonset_id'],$aint['geneset_id']);
		//print_r($aint);
		//preparing mysql connection and uploading info for each voucher
		foreach ($newsarr[0] as $vn){
			$prep = array();
			$qv = "SELECT * FROM temp_".$p_.$tablearraynames[$i]." WHERE ";
			if (count($tablearrayvitals[$i]) == 2) {
				$tav = explode("§",$vn);
				$qv .= $tablearrayvitals[$i][0]."='".$tav[0]."' AND " .$tablearrayvitals[$i][1]."='".$tav[1]."'";
			}
			else{
				$qv .= $tablearrayvitals[$i] . "='$vn'";
			}
			//echo "<br>qv=$qv<br>";print_r($tav[0]);echo"<br>";
			$rv = mysql_query($qv) or die("Error in query(7): $qv. " . mysql_error());
			//if records present
			if( mysql_num_rows($rv) > 0 ) {
				while( $rowv = mysql_fetch_object($rv) ) {
					foreach ($aint as $rva => $rvb){
						if ($rowv->$rvb != ''){
							$prep['values'][] = $rowv->$rvb;
							$prep['headers'][] = $rvb;
						}
					}
				}
			}
		//prepare and run INSERT statement
		$qvia = "INSERT INTO ". $p_ .$tablearraynames[$i]." (". implode(", ",$prep['headers']) .") VALUES ('". implode("', '",$prep['values']) ."')";
		$result = mysql_query($qvia) or die ("Error in query(8): $qvia. " . mysql_error());
		//insert mysql statement here
		if (count($tablearrayvitals[$i]) == 2) {
			$output0 .= "code: ".$tav[0]." for gene: ".$tav[1]."<br>";
		}
		else{
			$output0 .= "$vn<br>";
		}
		//echo "<br>$qvia<br>";
		}
	}
	// going through and updating records that exist but have additional info
	if (count($newsarr[1])>0){
		$output1 .= "<br>Updated ".$tablearraynames[$i]." (for former empty fields):<br>";
		foreach ($newsarr[1] as $vn1 => $vn1h){
			//echo "<br>$vn1<br>";print_r($vn1h);
			$hlist = implode(",",$vn1h);
			$qvn1h = array();
			$qvn1s = "SELECT $hlist FROM temp_" . $p_ .$tablearraynames[$i]." WHERE ";
			if (count($tablearrayvitals[$i]) == 2) {
				$tav = explode("§",$vn1);
				$qvn1s .= $tablearrayvitals[$i][0]."='".$tav[0]."' AND " .$tablearrayvitals[$i][1]."='".$tav[1]."'";
			}
			else{
				$qvn1s .= $tablearrayvitals[$i] . "='$vn1'";
			}
			$rvn1s = mysql_query($qvn1s) or die("Error in query(9): $qvn1s. " . mysql_error());
			//if records present
			if( mysql_num_rows($rvn1s) > 0 ) {
				while( $rowvn1s = mysql_fetch_object($rvn1s) ) {
					foreach ($vn1h as $h){
						$qvn1h[] = "$h='".$rowvn1s->$h."'";
					}
				}
			}
			$qvn1u = "UPDATE " . $p_ .$tablearraynames[$i]." SET ". implode(", ",$qvn1h) ." WHERE ";
			if (count($tablearrayvitals[$i]) == 2) {
				$tav = explode("§",$vn1);
				$qvn1u .= $tablearrayvitals[$i][0]."='".$tav[0]."' AND " .$tablearrayvitals[$i][1]."='".$tav[1]."'";
			}
			else{
				$qvn1u .= $tablearrayvitals[$i] . "='$vn1'";
			}
			// run query
			$result = mysql_query($qvn1u) or die ("Error in query(10): $qvn1u. " . mysql_error());
			
			if (count($tablearrayvitals[$i]) == 2) {
				$output1 .= "code: ".$tav[0]." for gene: ".$tav[1].": ". str_replace(",",", ",$hlist) ."<br>";
			}
			else{
				$output1 .= "$vn1: ". str_replace(",",", ",$hlist) ."<br>";
			}
			//echo "<br>$qvn1u<br>";
		}
	}
	// going through the records that are being left out due to conflict
	if (count($newsarr[2])>0){
		$output2 .= "<br> These ".$tablearraynames[$i]." will not be updated since they have conflicting values (printed) with already existing ".$tablearraynames[$i].":<br>";
		foreach ($newsarr[2] as $vn2 => $vn2h){
			if (count($tablearrayvitals[$i]) == 2) {
				$tav = explode("§",$vn2);
				$output2 .= "code: ".$tav[0]." for gene: ".$tav[1]." and values: ".implode(", ",$vn2h)."<br>";
			}
			else{
				if ($tablearraynames[$i] == 'genes' || $tablearraynames[$i] == 'vouchers'){
					$output2 .= "$vn2: ";
					foreach ($vn2h as $bby){
						if ($tablearraynames[$i] == 'genes'){$which_field = 'geneCode';}
						if ($tablearraynames[$i] == 'vouchers'){$which_field = 'code';}
						$q2 = "SELECT $bby FROM ".$p_.$tablearraynames[$i]." WHERE $which_field='$vn2'";
						$q2t = "SELECT $bby FROM temp_".$p_.$tablearraynames[$i]." WHERE $which_field='$vn2'";
						$r2 = mysql_query($q2) or die("Error in query(2b): $q2. " . mysql_error());
						$r2t = mysql_query($q2t) or die("Error in query(3): $q2t. " . mysql_error());
						$output2a =array();
						foreach (array($r2,$r2t) as $rnumb){
							if( mysql_num_rows($rnumb) > 0 ) {
							//if records present
								while( $row = mysql_fetch_object($rnumb) ) {
									$output2a[] .= $row->$bby;
								}
							}
							else{
								$output2a[] = "";
							}
						}
						$output2 .= "$bby: ('" . implode("' / '", $output2a) . "')  ";
					}
					$output2 .= "<br>";
				}
				else{
					$output2 .= "$vn2: " . implode(", ",$vn2h) ."<br>";
				}
			}
			$bad_vouchers[] = $vn2;
		}
	}
	if ($newsarr[3]>0){
		$output3 .= "<br> ". $newsarr[3] . " records from ".$tablearraynames[$i]." will not be updated since they already exist!<br>";
	}
}
if ($output_bads != '<br><b>These sequences/primers are not uploaded due to problems with connected codes/genes:</b><br>' || $output2 != '<br><b>These records are not updated due to conflict with existing records:</b><br>'){
	$output .= '<br><b>!!! There is some problems with the merging! See below under BAD!!!!</b><br><br><b>BAD!</b><br>';
	if ($output2 != '<br><b>These records are not updated due to conflict with existing records:</b><br>'){
		$output .= $output2;
	}
	if ($output_bads != '<br><b>These sequences/primers are not uploaded due to problems with connected codes/genes:</b><br>') {
		$output .= $output_bads;
	}
}
else{$output .= '<br><b>!!! Everything went ok with the merging!!!!</b><br>';}
$output .= "<br><b>GOOD!</b><br>".$output0 . $output1 . $output3;
//$output = str_replace(array("<br>","<b>","</b>","<u>","</u>"),array("\n","","","",""),$output);
echo "$output";
$tablearraynames = array("members", "vouchers", "genes", "sequences", "primers", "taxonsets", "genesets","search","search_results");
foreach ($tablearraynames as $taname){
	$q_droptemps = "DROP TABLE IF EXISTS `temp_".$p_.$taname."`";
	$r_droptemps = mysql_query($q_droptemps) or die ("Error in query(11): $q_droptemps. " . mysql_error());
}
if( file_exists($filename) ) {
	unlink($filename);
	echo "Finished importing database!";
}

?>