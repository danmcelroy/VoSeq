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


include'functions.php';
$admin = false;
$in_includes = false;
// need dojo?
$dojo = true;
// which dojo?
$whichDojo[] = 'Tooltip';

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
						Outgroup (code, for NEXUS and TNT): <input type="text" name="outgroup" size="10"><br />
						<img src="images/warning.png" /> Warning! your dataset will not necessarily be properly aligned! <img width="15px" height="16px" src="images/question.png" id="warn" alt="" /> 
								 <span dojoType="tooltip" connectId="warn" delay="1" toggle="explode">You need to be sure that your data is aligned!</span>
					</td>
				</tr>
				<tr>
					<td class="label"><b>Only for Protein coding:</br>(others will be output as standard)</b></br>Choose codon positions to use<br /><br /> (Override priority: <br />Amino Acids->Special->Degen->All->1st&#151;2nd,3rd):</td>
					<td class="field">
						Positions<br />
						&nbsp;&nbsp;&nbsp;<input type="checkbox" name="positions[all]" checked>all 
						&nbsp;&nbsp;&nbsp;<input type="checkbox" name="positions[1st]">1st 
						&nbsp;&nbsp;&nbsp;<input type="checkbox" name="positions[2nd]">2nd 
						&nbsp;&nbsp;&nbsp;<input type="checkbox" name="positions[3rd]">3rd<br /> 
						Partition by (positions)<br />
						&nbsp;&nbsp;&nbsp;<input type="radio" name="by_positions" value="asone" checked>as one 
						&nbsp;&nbsp;&nbsp;<input type="radio" name="by_positions" value="each">each  
						&nbsp;&nbsp;&nbsp;<input type="radio" name="by_positions" value="123">1st&#151;2nd, 3rd<br />
						Translations<br />
						&nbsp;&nbsp;&nbsp;<input type="checkbox" name="degen[yes]">Degen(erated) <img width="15px" height="16px" src="images/question.png" id="Degen" alt="" /> 
								 <span dojoType="tooltip" connectId="Degen" delay="1" toggle="explode">Degenerating nucleotides to IUPAC ambiguity codes at all sites that can potentially undergo<br />
											synonymous change in any and all pairwise comparisons of sequences in the data matrix,<br />
											thereby making synonymous change largely invisible and reducing the effect of<br />
											compositional heterogeneity but leaving the inference of non-synonymous change largely intact.<br /><br />
										From: "degen_v1_4.pl"; author: A. Zwick & A. Hussey; www.phylotools.com)<br /><br />Cite:<br />
										Zwick, A., Regier, J.C. & Zwickl, D.J. (2012). "Resolving Discrepancy between Nucleotides and<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
											Amino Acids in Deep-Level Arthropod Phylogenomics:Differentiating Serine Codons in 21-<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
											Amino-Acid Models". PLoS ONE 7(11): e47450.<br />
										Regier, J.C., Shultz, J.W., Zwick, A., Hussey, A., Ball, B., Wetzer, R. Martin, J.W. &<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
											Cunningham, C.W. (2010). "Arthropod relationships revealed by phylogenomic analysis of<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
											nuclear protein-coding sequences". Nature 463: 1079-1083.
								</span> &nbsp;&nbsp;=> as:
						&nbsp;&nbsp;<input type="radio" name="degenSZ" value="no" checked>normal 
						&nbsp;&nbsp;&nbsp;<input type="radio" name="degenSZ" value="S">S  <img width="15px" height="16px" src="images/question.png" id="S" alt="" /> 
								 <span dojoType="tooltip" connectId="S" delay="1" toggle="explode">Degen S: <br />alternative degen encoding of Serin 1 to Serin 2</span>
						&nbsp;&nbsp;&nbsp;<input type="radio" name="degenSZ" value="Z">Z <img width="15px" height="16px" src="images/question.png" id="Z" alt="" /> 
								 <span dojoType="tooltip" connectId="Z" delay="1" toggle="explode">Degen Z: <br />alternative degen encoding of Serin 2 to Serin 1</span>
						&nbsp;&nbsp;&nbsp;<input type="radio" name="degenSZ" value="SZ">SZ <img width="15px" height="16px" src="images/question.png" id="SZ" alt="" /> 
								 <span dojoType="tooltip" connectId="SZ" delay="1" toggle="explode">Degen SZ: <br />alternative degen encoding of Serin 1 and Serin 2 to NNN</span><br />
						&nbsp;&nbsp;&nbsp;<input type="checkbox" name="positions[aas]">Amino acids<br />
						Special<br />
						&nbsp;&nbsp;&nbsp;<input type="checkbox" name="positions[special]">Special <img width="15px" height="16px" src="images/question.png" id="Spec" alt="" /> 
								 <span dojoType="tooltip" connectId="Spec" delay="1" toggle="explode">Will enable you to specify different translations, positions and partitions<br />
									for each gene/alignment</span><br />
						
						
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
					<td class="label4" rowspan="2">
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
					<td class="label4" rowspan="2">
						Check to select</br> your alignment/gene:</br>
						</br>(if geneset chosen,</br> add extra genes)
					</td>
					<td class="field1" rowspan="2">
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
					
						
					<input type="submit" name="process_dataset" value="Create dataset" /><td>
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
