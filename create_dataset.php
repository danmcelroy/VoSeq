<?php
// #################################################################################
// #################################################################################
// Voseq create_dataset.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Script for input to dataset creation
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
#include('functions.php');

// print beginning of html page -- headers
include_once'includes/header.php';
nav();
// #################################################################################
// Section: Quering DB for genecodes and taxon sets
// #################################################################################
// open database connections
@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect');
mysql_select_db($db) or die ('Unable to select database');
$result = mysql_query("set names utf8") or die("Error in query: $query. " . mysql_error());
$gCquery = "SELECT geneCode FROM ". $p_ . "genes ORDER BY geneCode";
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
// #################################################################################
// Section: Output
// #################################################################################
// begin HTML page content
echo "<div id=\"content\">";
?>

<b>Select sequences you want for a dataset by entering the voucher codes and gene codes.</b> <a href="http://carlosp420.github.io/VoSeq/#create-datasets">See more help.</a>

<form action="includes/process_dataset.php" method="post">
<table border="0" width="960px" cellpadding="5px"> <!-- super table -->
	<tr>
		<td valign="top" colspan="2">
			<h1>Create dataset</h1>
			<table border="0" width="800px" cellspacing="0" cellpadding="0">
				<caption>Enter the required info to make yourself a ready-to-run dataset</caption>
				<tr>
					<td class="label">Choose file format:</td>
					<td class="field">
						<input type="radio" name="format" value="TNT">TNT format<br />
						<input type="radio" name="format" value="NEXUS">NEXUS format<br />
						<input type="radio" name="format" value="PHYLIP">PHYLIP format<br />
						<input type="radio" name="format" value="FASTA" checked>Unaligned FASTA format
						&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
						<br />
						Outgroup (code, for NEXUS and TNT): <input type="text" name="outgroup" size="10">
					</td>
				</tr>
				<tr>
					<td class="label"><b>Only for Protein coding:</br>(others will be output as standard)</b></br>Choose codon positions to use<br /><br /> (Override priority: <br />Amino Acids->Special->All->1st&#151;2nd,3rd):</td>
					<td class="field">
						<input type="checkbox" name="positions[all]" checked>all 
						<input type="checkbox" name="positions[1st]">1st 
						<input type="checkbox" name="positions[2nd]">2nd 
						<input type="checkbox" name="positions[3rd]">3rd 
						<input type="checkbox" name="positions[special]">Special
						<input type="checkbox" name="positions[aas]">Amino acids<br />
						Partition by (positions):
						<input type="radio" name="by_positions" value="asone" checked>as one 
						<input type="radio" name="by_positions" value="each">each  
						<input type="radio" name="by_positions" value="123">1st&#151;2nd, 3rd<br />
						<img src="images/warning.png" /> Warning! your dataset will not necessarily be properly aligned!
					</td>
				</tr>
				<tr>
					<td class="label">
						What info do you want in the taxon names?
					</td>
					<td class="field">
						<table>
						<td><input type="checkbox" name="taxonadds[code]"checked>Code</td>
						<td><input type="checkbox" name="taxonadds[orden]">Order</td>
						<td><input type="checkbox" name="taxonadds[family]">Family</td>
						<td><input type="checkbox" name="taxonadds[subfamily]">Subfamily</td>
						<td><input type="checkbox" name="taxonadds[tribe]">Tribe</td>
						<td><input type="checkbox" name="taxonadds[subtribe]">Subtribe</td>

						</tr />

						<td><input type="checkbox" name="taxonadds[genus]"checked>Genus</td>
						<td><input type="checkbox" name="taxonadds[species]"checked>Species</td>
						<td><input type="checkbox" name="taxonadds[subspecies]">Subspecies</td>
						<td><input type="checkbox" name="taxonadds[auctor]">Auctor</td>
						<td><input type="checkbox" name="taxonadds[hostorg]">Host org.</td>
						<td><input type="checkbox" name="taxonadds[genecode]">Gene code</td>
						</table>
					</td>
				</tr>
				<tr>
					<td class="label">
						For single gene datasets, </br>exclude taxa missing this gene?
					</td>
					<td class="field">
						<table>
						<td><input type="radio" name="exclude_missing" value="yes" checked>yes</td>
						<td><input type="radio" name="exclude_missing" value="no" >no</td>
						</table>
					</td>
				</tr>
				<tr>
					<td class="label">
						For multigene datasets, </br>exclude taxa with less than X genes?
					</td>
					<td class="field">
						<table>
							Minimum number of genes: <input type="text" name="less_than_genes" size="3">
						</table>
					</td>
				</tr>
				<tr>
					<td class="label">
						Ignore introns? </br>('yes' will not use them in the data set)
					</td>
					<td class="field">
						<table>
						<td><input type="radio" name="ignore_introns" value="yes" checked>yes</td>
						<td><input type="radio" name="ignore_introns" value="no" >no</td>
						</table>
					</td>
				</tr>
			</table>
		</td>
	</tr>
	<tr>
		<td align="left" colspan="2">
			<table border="0" cellspacing="5px">
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
					</td>
				</tr>
				<tr>
					<td class="label4">
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
						&nbsp;&nbsp;&nbsp;&nbsp;and so on...
					</td>
					<td>
					</td>
					<td class="label4">
						Check to select</br> your alignment/gene:
					</td>
					<td class="field1">
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
					</td>
					<td>
						<textarea name="codes" rows="10">
						</textarea>
					</td>
					<td></td><td><input type="submit" name="process_dataset" value="Create dataset" /><td>
				</tr>
				<tr>
					<td>
					</td>
				</tr>
			</table>
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
