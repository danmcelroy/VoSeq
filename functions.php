<?php
// #################################################################################
// #################################################################################
// Voseq functions.php
// author(s): Carlos Pe�a & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Certain funtions used as include() in other scripts
// #################################################################################
// #################################################################################
// Section: Startup/includes
// #################################################################################
// includes
ob_start();//Hook output buffer - disallows web printing of file info...
include'conf.php';
ob_end_clean();//Clear output buffer//includes

// #################################################################################
// Section: array_combine() function for PHP versions earlier than 5 
// #################################################################################
if (!function_exists('array_combine')) // for php 4 
{
	function array_combine($arr1,$arr2) {
		$out = array();
		foreach ($arr1 as $key1 => $value1) {
			$out[$value1] = $arr2[$key1];
		}
		return $out;
	}
}

// #################################################################################
// Section: stripos() function for PHP versions earlier than 5 
// #################################################################################
if (!function_exists('stripos')) // for php 4
{
	function stripos($haystack, $needle){
		return strpos($haystack, stristr( $haystack, $needle ));
	}
}


// #################################################################################
// Section: story.php Shows voucher images in Voucher page (multi photos)
// It will show up to two photos on top right
// All other photos will be shown at the bottom with another function
// #################################################################################
function show_multi_photos($flickr_id, $voucherImage, $thumbnail, $admin) {
	if( $voucherImage == NULL && $thumbnail == NULL ) {
		echo "<div class='voucher'>";
		echo "<img class=\"voucher\" src=\"na.gif\"/>";
		echo "</div>";
	}
	else {
		$f = explode("|", $flickr_id);
		$v = explode("|", $voucherImage);
		$t = explode("|", $thumbnail);

		# create an array with each photo data as array as well 
		$j = 0;
		$photos = array();
		while( $j < count($v) ) {
			if( $j < 3 ) {
				$tmp = new stdClass;
				$tmp->flickr_id = trim($f[$j]);
				$tmp->voucherImage = trim($v[$j]);
				$tmp->thumbnail = trim($t[$j]);

				$photos[] = $tmp;
			}
			$j = $j + 1;
		}

		$i = 1;
		foreach( $photos as $photo ) {
			$photo->flickr_id = str_replace("|", "", $photo->flickr_id);
			$photo->flickr_id = trim($photo->flickr_id);

			if( $photo->voucherImage != "" && $photo->thumbnail != "" ) {
				if( $admin == true ) {
					echo "\n<div class=\"voucher\" id=\"" . $i . "\">";
					echo "\n<a href='#' title='Delete photo' class='delete'>";
					echo "<img class='delete' src='images/delete.png' /></a>";
					echo "\n<a href='" . $photo->voucherImage . "' target=\"_blank\">";
					echo "\n<img class='voucher' src=\"" . $photo->thumbnail . "\" />";
					echo "\n</a>";
					echo "\n</div>";
					echo "<br />";
				}
				else {
					echo "\n<div class=\"voucher\">";
					echo "\n<a href='" . $photo->voucherImage . "' target='_blank'>";
					echo "\n<img class='voucher' src=\"" . $photo->thumbnail . "\"/>";
					echo "\n</a>";

					if( $photo->flickr_id != "" && $photo->flickr_id != NULL ) {
						echo "\n<div class='eol_button' onclick='send_to_EOL(". $photo->flickr_id . ");'>";
						echo "<img src='images/eol_button.png' ";
						echo "id='" . $photo->flickr_id . "' alt='' />";
						echo "Share photo with EOL</div>";
					}

					echo "</div>\n";
					echo "<br />";
				}
			}
			$i = $i + 1;
		}
	}
}

function show_all_other_photos($flickr_id, $voucherImage, $thumbnail, $admin) {
	if( $voucherImage == NULL && $thumbnail == NULL ) {
		echo "";
	}
	elseif ( $voucherImage == "na.gif" || $thumbnail == "na.gif" ) {
		echo "";
	}
	else {
		$f = explode("|", $flickr_id);
		$v = explode("|", $voucherImage);
		$t = explode("|", $thumbnail);

		# show from photo 3 onwards
		if( count($v) > 3 ) {
			$j = 3;
			$photos = array();
			while( $j < count($v) ) {
				if( $j > 2 ) {
					$tmp = new stdClass;
					$tmp->flickr_id = trim($f[$j]);
					$tmp->voucherImage = trim($v[$j]);
					$tmp->thumbnail = trim($t[$j]);

					$photos[] = $tmp;
				}
				$j = $j + 1;
			}
			
			echo "<tr>\n<td>\n";
			echo "<table border='0' width='760px' cellpadding='5px'>";
			echo "<th><h1>Additional photos</h1></th>";
			echo "<tr>";
			$i = 3;
			foreach( $photos as $photo ) {
				$photo->flickr_id = str_replace("|", "", $photo->flickr_id);
				$photo->flickr_id = trim($photo->flickr_id);
	
				if( $photo->voucherImage != "" && $photo->thumbnail != "" ) {
					echo "<td width='200px'>\n";
					echo "<div class='voucher' id=\"". $i . "\">\n";

					if( $admin == true ) {
						echo "<a href='#' class='delete'><img class='delete'";
						echo " src='images/delete.png' title='Delete photo' /></a>\n";
					}
					echo "<a href='" . $photo->voucherImage . "' target=\"_blank\">";
					echo "<img class=\"voucher\" src=\"" . $photo->thumbnail . "\"/></a>\n";

					if( $photo->flickr_id != "" && $photo->flickr_id != NULL ) {
						echo "\n<div class='eol_button' onclick='send_to_EOL(". $photo->flickr_id . ");'>";
						echo "<img src='images/eol_button.png' ";
						echo "id='" . $photo->flickr_id . "' alt='' />";
						echo "Share photo with EOL</div>";
					}
					echo "</div>\n";
					echo "</td>\n";
				}
				if( $i % 3 == 2 ) {
					echo "</tr><tr>";
				}
				$i = $i + 1;
			}
			echo "\n</tr>\n";
			echo "</table>";
			echo "</td>\n</tr>\n";
		}
	}
}



// #################################################################################
// Section: dofastafiles() function 
// returns a fasta file in GenBank format
// #################################################################################
function dofastafiles($geneCode, $code, $p_) {
	// open file
	$cwd = getcwd();
	$fastafile = $cwd . '/myfastafile.txt';
	$genbank_fastafile = $cwd . '/my_genbank_fastafile.txt';
	$handle1 = fopen($fastafile, "a");
	$handle2 = fopen($genbank_fastafile, "a");
	$query3 = "SELECT " . $p_ . "vouchers.orden, 
					  " . $p_ . "vouchers.family, 
					  " . $p_ . "vouchers.subfamily, 
					  " . $p_ . "vouchers.tribe, 
					  " . $p_ . "vouchers.genus, 
					  " . $p_ . "vouchers.species, 
					  " . $p_ . "sequences.sequences, 
					  " . $p_ . "genes.description FROM 
					  " . $p_ . "vouchers, 
					  " . $p_ . "sequences, 
					  " . $p_ . "genes WHERE 
					  " . $p_ . "vouchers.code='$code' AND 
					  " . $p_ . "sequences.code='$code' AND 
					  " . $p_ . "sequences.geneCode='$geneCode' AND 
					  " . $p_ . "genes.geneCode='$geneCode'";
	$result3 = mysql_query($query3) or die("Error in query: $query3. " . mysql_error());
	if (mysql_num_rows($result3) > 0) {
		$lineage = " [Lineage=Eukaryota; Metazoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota;";
		while ($row3 = mysql_fetch_object($result3)) {
			if( $row3->orden ) {
				$lineage .= " $row3->orden;";
			}
			if( $row3->family ) {
				$lineage .= " $row3->family;";
			}
			if( $row3->subfamily ) {
				$lineage .= " $row3->subfamily;";
			}
			if( $row3->tribe ) {
				$lineage .= " $row3->tribe;";
			}
			$lineage .= " $row3->genus] ";
			$species = str_ireplace(" ", "_", $row3->species);
			fwrite($handle1, ">" . $row3->genus . "_" . $species . "_" . $code . "\n$row3->sequences\n");
			fwrite($handle2, ">" . $row3->genus . "_" . $species . "_" . $code . " [org=$row3->genus $row3->species] [Specimen-voucher=$code]");
			fwrite($handle2, " [note=" . $row3->description . "] $lineage");
			fwrite($handle2, "\n$row3->sequences\n");
		}
		unset($lineage);
	}
	fclose($handle1);
	fclose($handle2);
}
// #################################################################################
// Section: dofastafile() function 
// returns a fasta file in standard format
// #################################################################################
function dofastafile($geneCode, $code, $p_){
					// 			open file
$cwd = getcwd();
$fastafile = $cwd . '/myfastafile.txt';
$handle = fopen($fastafile, "a");
$query3 = "SELECT " . $p_ . "vouchers.code, 
				  " . $p_ . "vouchers.genus, 
				  " . $p_ . "vouchers.species, 
				  " . $p_ . "sequences.sequences FROM 
				  " . $p_ . "vouchers, 
				  " . $p_ . "sequences WHERE 
				  " . $p_ . "vouchers.code='$code' AND 
				  " . $p_ . "sequences.code='$code' AND geneCode='$geneCode'";
$result3 = mysql_query($query3) or die("Error in query: $query3. " . mysql_error());
if (mysql_num_rows($result3) > 0)
	{
	while ($row3 = mysql_fetch_object($result3))
		{
		fwrite($handle, ">" . $row3->genus . "_" . $row3->species . "_" . $row3->code . "\n$row3->sequences\n");
		}
	}
fclose($handle);
}
// #################################################################################
// Section: formatdate() function 
// format MySQL DATETIME value into a more readable string
// #################################################################################
// 
function formatDate($val, $date_timezone, $php_version) {
	if( $php_version == "5" ) {
		date_default_timezone_set($date_timezone); //php5
	}
	$arr = explode('-', $val);
	$day = explode(" ", $arr[2]);
	return date('d M Y', mktime(0,0,0, $arr[1], $day[0], $arr[0]));
}
// #################################################################################
// Section: getSeqs() function 
// get and print all sequences for a voucher, for the "story" (story.php) page
// #################################################################################
function getSeqs($code, $host, $user, $pass, $db, $p_) {
	// open database connection
	$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect! ' . mysql_error());
	// select database
	mysql_select_db($db) or die ('Unable to select database ' . mysql_error());
	if( function_exists(mysql_set_charset) ) {
		mysql_set_charset("utf8");
	}
	// generate and execute query
	$query  = "SELECT geneCode, CHAR_LENGTH(sequences) AS mylen, (2*CHAR_LENGTH(sequences) - CHAR_LENGTH(REPLACE(sequences, '?', '')) - CHAR_LENGTH(REPLACE(sequences, '-', ''))) AS amb, LEFT((labPerson),8) AS labPerson, accession FROM " . $p_ . "sequences WHERE code='$code' ORDER BY geneCode";
	
	$result  = mysql_query($query)  or die("Error in query: $query.  " . mysql_error());
	ob_start();//Hook output buffer - disallows web printing of file info...
	include 'conf.php';
	ob_end_clean();//Clear output buffer//includes
	if (mysql_num_rows($result) > 0) {
		$i = "0";
		while ($row = mysql_fetch_object($result)) {
			$i += "1";
			$labPerson = explode(" ", $row->labPerson);
			$labPerson = $labPerson[0];
			echo "<tr>";
			// masking URLs, this variable is set to "true" or "false" in conf.php file
			if($mask_url =="true") {
				echo "<td class=\"field4\"><a href='" . $base_url . "/home.php'  onclick=\"return redirect('sequences.php?code=". $code . "&amp;geneCode=" . $row->geneCode . "')\">" . $row->geneCode . "</a></td>";
			}
			else {
				echo "<td class=\"field4\"><a href=\"sequences.php?code=". $code . "&amp;geneCode=" . $row->geneCode . "\">" . $row->geneCode . "</a></td>";
			}
			echo "<td class=\"field4\">" . $row->mylen . "</td>";
			echo "<td class=\"field4\">" . $row->amb . "</td>";
			echo "<td class=\"field4\">" . $labPerson . "</td>";
			echo "<td class=\"field4\"><a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=nucleotide&amp;term=" . $row->accession . "[accn]&amp;doptcmdl=GenBank\" target=\"_blank\">" . $row->accession . "</a>&nbsp;</td>";

			echo "<td class=\"field4\">";
			echo "<a href=\"includes/blast_locally.php?code=". $code . "&amp;geneCode=$row->geneCode\" target=\"_blank\" ><img class=\"link\" width=\"18px\" height=\"14px\" src=\"images/local_blast.png\" id=\"local_blast" . $i . "\"/></a>";	
			echo "<span dojoType=\"tooltip\" connectId=\"local_blast" . $i . "\" delay=\"1\" toggle=\"explode\">Local blast: this sequence against <br /><b>all sequences of the same gene code</b></span>";

			echo "&nbsp;<a href=\"includes/blast_locally_full_db.php?code=". $code . "&amp;geneCode=$row->geneCode\" target=\"_blank\" ><img class=\"link\" width=\"16px\" height=\"16px\" src=\"images/database.png\" id=\"full_database" . $i . "\" /></a>";	
			echo "<span dojoType=\"tooltip\" connectId=\"full_database" . $i . "\" delay=\"1\" toggle=\"explode\"><b>FULL blast:</b> this sequence against <b>the full database</b></span>";
			echo "</td>";

			echo "<td class=\"field5\">";
			$tmp = explode("-", $row->geneCode);
			echo "<a href=\"includes/blast_vs_genbank.php?code=". $code ."&amp;geneCode=$row->geneCode\" target=\"_blank\"><img class=\"link\" width=\"18px\" height=\"14px\" src=\"images/ncbi_blast.png\" id=\"vs_genbank" . $i . "\" /></a>";
			echo "<span dojoType=\"tooltip\" connectId=\"vs_genbank" . $i . "\" delay=\"1\" toggle=\"explode\">BLAST against GenBank</span>";
			unset($tmp);
			echo "</td>\n\t\t</tr>\n\t\t";
  		}
		// print empty fields and finish this table
		echo "<tr>\n\t\t\t";
		//link to add new sequence for this record
		echo "<td class=\"field4\"><a href=\"admin/listseq.php?code=" . $code . "\">add seq</a></td>";
		echo "<td class=\"field4\">&nbsp;</td>
				<td class=\"field4\">&nbsp;</td>
				<td class=\"field4\">&nbsp;</td>
				<td class=\"field4\">&nbsp;</td>
				<td class=\"field4\">&nbsp;</td>
				<td class=\"field5\">&nbsp;</td>";
		echo "</tr></table>\n\t<!-- end colmun 2 second block -->\n";
		}
	else
		{
		echo "<tr>
			<td class=\"field4\"><a href=\"admin/listseq.php?code=" . $code . "\">Add Seq</a></td>
			<td class=\"field4\">&nbsp;</td>
			<td class=\"field4\">&nbsp;</td>
			<td class=\"field4\">&nbsp;</td>
			<td class=\"field4\">&nbsp;</td>
			<td class=\"field4\">&nbsp;</td>
			<td class=\"field5\">&nbsp;</td>
			</tr>
			</table><!-- 	end column 2, second block -->\n\n";
		}
	}

// #################################################################################
// Section: clean_string() function 
// returns a string cleaned of certain characters
// #################################################################################
function clean_string($string) {
	$i = 0;
	if( (isset($string) && trim($string) != '') ) {
		$user_strings = array();
		$symbols = array(",",'"',"'","&","/","\\",";","=");
		#is number? then dont filter by symbols
		if( is_numeric($string) ) {
			array_push($user_strings, $string);
		}
		else { #not number, then clean by filtering symbols
			$id_subject = trim(str_replace($symbols, "", $string));
			$subject = explode(" ", $id_subject);
			foreach( $subject as $val ) {
				if( trim($val) != "" ) {
					$pattern = '/[a-�A-�0-9_\.\-]+/';
					preg_match($pattern, $val, $match);
					if( $i < 3 ) {
						array_push($user_strings, $match[0]);
					}
					$i++;
				}
			}
		}
		return $user_strings;
	}
}

// #################################################################################
// Section: update_comboboxes() function 
// updates the Dojo comboboxes in /Dojo
// #################################################################################
function update_comboboxes() {

	ob_start();//Hook output buffer - disallows web printing of file info...
	include 'conf.php';
	ob_end_clean();//Clear output buffer//includes

	/*** make dojos .js files ***/

	// Initialize default settings
	$MYSQL_PATH = '/usr/local/mysql/bin';

    $comboName[] = 'orden';
	$comboName[] = 'family';
	$comboName[] = 'subfamily';
	$comboName[] = 'tribe';
    $comboName[] = 'subtribe';
	$comboName[] = 'genus';
    $comboName[] = 'species';
	$comboName[] = 'subspecies';
	$comboName[] = 'country';
    $comboName[] = 'code';
	$comboName[] = 'collector';
	$comboName[] = 'extractor';
	$comboName[] = 'hostorg';

	// table sequences
    $comboNameSeq[] = 'labPerson';
	$comboNameSeq[] = 'geneCode';
	$comboNameSeq[] = 'accession';

    // table primers
	$comboNamePri[] = 'primer1';
	$comboNamePri[] = 'primer2';
	$comboNamePri[] = 'primer3';
    $comboNamePri[] = 'primer4';
	$comboNamePri[] = 'primer5';
    $comboNamePri[] = 'primer6';
	
	$cwd = dirname(__FILE__);

    // connect to database
	@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect ' . mysql_error());
	mysql_select_db($db) or die ('Unable to select database '. mysql_error());
	if( function_exists(mysql_set_charset) ) {
    	mysql_set_charset("utf8");
	}

	// do table vouchers
	foreach ($comboName as $value) {
		if( !is_dir($local_folder . '/dojo_data') ) {
			mkdir($local_folder . '/dojo_data', 0755);
		}

    	$query = "SELECT DISTINCT $value FROM " . $p_ . "vouchers ORDER BY $value ASC";
		$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
		
		$comboFile = $local_folder . '/dojo_data/comboBoxData_' . $value . '.js';
    	if ( file_exists($comboFile) ) {
    		unlink($comboFile);
		}
		$handle = fopen($comboFile, "w");
    	fwrite($handle, "[\n");
		
		while( $row = mysql_fetch_object($result) ) {
    		if ( $row->$value == "" ) {
    			continue;
			}
			else {
				fwrite($handle, "\t[\"" . $row->$value . "\"],\n");
			}
		}
    	fwrite($handle, "]\n");
		//echo "$value\n";
		fclose($handle);
	}

	// do table sequences
    foreach ($comboNameSeq as $value) {
		$query = "SELECT DISTINCT $value FROM " . $p_ . "sequences ORDER BY $value ASC";
    	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
		
		$comboFile = $local_folder . '/dojo_data/comboBoxData_' . $value . '.js';
		if ( file_exists($comboFile) ) {
			unlink($comboFile);
    	}
		$handle = fopen($comboFile, "w");
		fwrite($handle, "[\n");
	
    	while( $row = mysql_fetch_object($result) ) {
			if ( $row->$value == "" ) {
    			continue;
			}
    		else {
				fwrite($handle, "\t[\"" . $row->$value . "\"],\n");
			}
		}
    	fwrite($handle, "]\n");
		//echo "$value\n";
		fclose($handle);
	}

	// do table primers
    foreach ($comboNamePri as $value) {
		$query = "SELECT DISTINCT $value FROM ". $p_ . "primers ORDER BY $value ASC";
    	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
		
		$comboFile = $local_folder . '/dojo_data/comboBoxData_' . $value . '.js';
		if ( file_exists($comboFile) ) {
			unlink($comboFile);
    	}
		$handle = fopen($comboFile, "w");
		fwrite($handle, "[\n");
	
    	while( $row = mysql_fetch_object($result) ) {
			if ( $row->$value == "" ) {
    			continue;
			}
    		else {
    				fwrite($handle, "\t[\"" . $row->$value . "\"],\n");
			}
		}
		fwrite($handle, "]\n");
    	//echo "$value\n";
		fclose($handle);
    }
}


// #################################################################################
// Section: includes/process_dataset.php
// has_seqs? checks if array of sequences consist only of "?"
// @input: $seqout_array
// @output: string "true" or "false"
// #################################################################################
function has_seqs($seqout_array, $gen) {
	$total_seqs = "";

	if( count($seqout_array) > 0 ) {
		foreach($seqout_array[$gen] as $code => $seq) {
			$total_seqs .= $seq;
		}
	}
	$total_seqs = str_ireplace("?", "", $total_seqs);
	$total_seqs = str_ireplace("N", "", $total_seqs);
	$total_seqs = str_ireplace("-", "", $total_seqs);
	$total_seqs = trim($total_seqs);
	if( strlen($total_seqs) > 0) {
		return "true";
	}
	else {
		return "false";
	}
}


// #################################################################################
// Section: 
// does a curl request
// @input: $url string
// @output: result using return transfer as string
// #################################################################################
function get_from_URL($url){
	// use cURL, safer than simplexml_load_file
	$ch = curl_init();

	curl_setopt($ch, CURLOPT_URL, $url);
	curl_setopt($ch, CURLOPT_HEADER, 0);
	curl_setopt($ch, CURLOPT_RETURNTRANSFER, TRUE);
	curl_setopt ($ch, CURLOPT_USERAGENT, 'Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 5.1; .NET CLR 1.1.4322)');
	curl_setopt ($ch, CURLOPT_SSL_VERIFYPEER, false);
	curl_setopt ($ch, CURLOPT_SSL_VERIFYHOST, false);

	$result = curl_exec($ch);

	curl_close($ch);

	return $result;
}


// #################################################################################
// Section: admin
// cleans fields voucherImage, thumbnail for photos, removes double | and makes sure
// there is always one | at the begining of field
// @input: $p_ table prefix
// @input: voucher id
// #################################################################################
function clean_fields($p_, $id) {
	$query2  = "SELECT voucherImage, flickr_id, thumbnail FROM ". $p_;
	$query2 .= "vouchers WHERE id='" . $id . "'";
	$result2 = mysql_query($query2) or die("Error in query: $query2 ". mysql_error());
	if( mysql_num_rows($result2) > 0 ) {
		while( $row2 = mysql_fetch_object($result2) ) {
			$v = trim($row2->voucherImage);
			$f = trim($row2->flickr_id);
			$t = trim($row2->thumbnail);

			$v = preg_replace("/^\|{2,}/", "|", $v);
			$f = preg_replace("/^\|{2,}/", "|", $f);
			$t = preg_replace("/^\|{2,}/", "|", $t);

			$v = preg_replace("/\|+/", "|", $v);
			$f = preg_replace("/\|+/", "|", $f);
			$t = preg_replace("/\|+/", "|", $t);

			$v = preg_replace("/^\|$/", "", $v);
			$f = preg_replace("/^\|$/", "", $f);
			$t = preg_replace("/^\|$/", "", $t);

			$q2  = "UPDATE ". $p_ . "vouchers SET ";
			$q2 .= "voucherImage='$v', flickr_id='$f', thumbnail='$t' ";
			$q2 .= "WHERE id=$id";	
			mysql_query($q2) or die("Error in query: $q2 ". mysql_error());
		}
	}
}

// #################################################################################
// Section: Remove introns
// removes introns from sequence string
// @input: $sequence and $introns from genetable (w-x;y-z...)
// @output: result shortened sequence, new length, number of introns, total intron 
//          length and intron sequences in array using return transfer
// #################################################################################
function remove_introns ($sequence, $intron_string){
	$introns_rev  = array_reverse(explode(';',$intron_string));
	$numtrons = count($introns_rev);
	$seq_array = preg_split('#(?<=.)(?=.)#s', $sequence); // making sequence/nucleotide array
	$intron_array;
	$intron_total_l = 0;
	foreach ($introns_rev as $inrev){
		$ir = explode('-',$inrev);
		$len = $ir[1] - $ir[0] + 1;
		$intron_total_l = $intron_total_l + $len;
		$start = $ir[0]-1;
		$intron_array[] = implode('',array_splice($seq_array, $start, $len));
	}
	$intron_array = array_reverse($intron_array);
	$seqlen = count($seq_array);
	$seqout = implode('',$seq_array);
	$array_out = array($seqout,$seqlen,$numtrons, $intron_total_l);
	foreach ($intron_array as $ia){$array_out[] = $ia;}
	return $array_out;
}

?>