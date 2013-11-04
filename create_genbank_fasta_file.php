<?php
// #################################################################################
// #################################################################################
// Voseq create_genbank_fasta_file.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Script for input to GeneBank FASTA file creation
// #################################################################################
// #################################################################################
// Section: Startup/includes
// #################################################################################
//check login session
include'login/auth.php';
// includes
include'markup-functions.php';
ob_start();//Hook output buffer - disallows web printing of file info...
include 'conf.php';
ob_end_clean();//Clear output buffer//includes

// print beginning of html page -- headers
include_once'includes/header.php';
nav();
// #################################################################################
// Section: Quering DB for genecodes and taxonsets
// #################################################################################
// open database connections
@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect');
mysql_select_db($db) or die ('Unable to select database');
if( function_exists(mysql_set_charset) ) {
	mysql_set_charset("utf8");
}

$gCquery = "SELECT geneCode FROM " . $p_ . "genes ORDER BY geneCode";
$gCresult = mysql_query($gCquery) or die("Error in query: $query. " . mysql_error());
// if records present
$geneCodes_array = array();
if( mysql_num_rows($gCresult) > 0 ) {
	while( $row = mysql_fetch_object($gCresult) ) {
		$geneCodes_array[] = $row->geneCode;
	}
}
// create dataset list
$Tsquery = "SELECT taxonset_name FROM ". $p_ . "taxonsets ORDER BY taxonset_name";
$Tsresult = mysql_query($Tsquery) or die ("Error in query: $Tsquery. " . mysql_error());
// if records present
$taxonsets = array();
if( mysql_num_rows($gCresult) > 0 ) {
	while( $rowTs = mysql_fetch_object($Tsresult) ) {
		$taxonsets[] = $rowTs->taxonset_name;
	}
}
else {unset($taxonsets);}

// create geneset list
$Tsquery = "SELECT geneset_name FROM ". $p_ . "genesets ORDER BY geneset_name";
$Tsresult = mysql_query($Tsquery) or die ("Error in query: $Tsquery. " . mysql_error());
// if records present
$genesets = array();
if( mysql_num_rows($gCresult) > 0 ) {
	while( $rowTs = mysql_fetch_object($Tsresult) ) {
		$genesets[] = $rowTs->geneset_name;
	}
}
else {unset($taxonsets);}
// #################################################################################
// Section: Output
// #################################################################################
// begin HTML page content
echo "<div id=\"content\">";
?>

<form action="includes/make_fasta_genbank.php" method="post">
<h1>Create GenBank FASTA file</h1>
<table border="0" width="800px" cellpadding="0px" cellspacing="5"> <!-- super table -->
	<caption>Enter the required info to make yourself a FASTA file to be submitted to GenBank:</caption>
	<tr>
		<td class="label4">
			Choose taxonset
		</td>
		<td class="field1">
			<select name="taxonsets" size="1" style=" BORDER-BOTTOM: outset; BORDER-LEFT: 
			outset; BORDER-RIGHT: outset; BORDER-TOP: outset; FONT-FAMILY: 
			Arial; FONT-SIZE: 12px"> 
			<option selected value="Choose taxonset">Choose taxonset</option> 
			<?php  // create a pulldown-list with all taxon set names in the db
			if (isset($taxonsets)){
			foreach ($taxonsets as $taxonset){ echo "<option value=\"$taxonset\">$taxonset</option> ";}
			}
			?>
			</select>
		</td><td></td>
		<td class="label4">
			Choose geneset
		</td>
		<td class="field1">
			<select name="genesets" size="1" style=" BORDER-BOTTOM: outset; BORDER-LEFT: 
			outset; BORDER-RIGHT: outset; BORDER-TOP: outset; FONT-FAMILY: 
			Arial; FONT-SIZE: 12px"> 
			<option selected value="Choose geneset">Choose geneset</option> 
			<?php  // create a pulldown-list with all geneset names in the db
			if (isset($genesets)){
			foreach ($genesets as $geneset){ echo "<option value=\"$geneset\">$geneset</option> ";}
			}
			?>
			</select>
		</td>
	</tr>
	<tr>
		<td class="label4" rowspan=2>
			...and/or a <br />list of codes:
		</td>
		<td class="field1">
			Enter one code per line</br>in the box below.</br>
			(with a -- sign before taxon</br>
			names to disable them from </br>the taxon set)</br>
			For example:<br />
			&nbsp;&nbsp;&nbsp;&nbsp;tA1<br />
			&nbsp;&nbsp;&nbsp;&nbsp;S077<br />
			&nbsp;&nbsp;&nbsp;&nbsp;--S078<br />
			&nbsp;&nbsp;&nbsp;&nbsp;and so on...<br />
		
		<textarea name="codes" rows="10">
			</textarea></td>
		<td>
		</td>
		<td class="label4" rowspan=2>
			Check to select</br> your alignment/gene:</br>
			</br>(if geneset chosen,</br> add extra genes)
		</td>
		<td class="field1" rowspan=2>
			<?php $i = 0;
					echo "<table>";
					foreach ($geneCodes_array as $genes) {
						$i = $i +1;
						echo "<td><input type=\"checkbox\" name=\"geneCodes[$genes]\" />$genes</td>"; 
						if ($i == 4) {
							echo "</tr />";
							$i = 0;
						}
					}
					echo "</table>";
			?>
		</td>
	</tr>
	<tr>
		<td>
		
			
		<input type="submit" name="make_table" value="Make genBank table" />
	</tr>
</table>

</form>

</div> <!-- end content -->

<?php
make_footer($date_timezone, $config_sitename, $version, $base_url);
?>
</body>
</html>
