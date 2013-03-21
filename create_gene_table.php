<?php
// #################################################################################
// #################################################################################
// Voseq create_table.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Input for creation of gene (xml) tables
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
#include('includes/yahoo_map.php');
#include('includes/show_coords.php');

// #################################################################################
// Section: Query DB for genecodes and taxonsets
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
//create dataset list
$Tsquery = "SELECT taxonset_name FROM " . $p_ . "taxonsets ORDER BY taxonset_name";
$Tsresult = mysql_query($Tsquery) or die ("Error in query: $Tsquery. " . mysql_error());
// if records present
$taxonsets = array();
if( mysql_num_rows($gCresult) > 0 ) {
	while( $rowTs = mysql_fetch_object($Tsresult) ) {
		$taxonsets[] = $rowTs->taxonset_name;
	}
}
else {unset($taxonsets);}

// print beginning of html page -- headers
include_once'includes/header.php';
nav();

// #################################################################################
// Section: Output
// #################################################################################
// begin HTML page content
echo "<div id=\"content\">";
?>

<b>You can create a MS Excel table of genes/alignments for a certain dataset or list of taxa. <br /> 
Instead of typing your specimen codes in the text area below, you could select a Taxonset  
(provided that it has been set <a href="admin/add_taxonset.php">here</a>).<br /></br>
This can only be done for aligned genes/alignments!</br></br>
This table will be ready to attach to a manuscript for publication.</b>


<form action="includes/make_gene_table.php" method="post">
<table border="0" width="960px" cellpadding="5px"> <!-- super table --> 
	<tr>
		<td valign="top" colspan="2">
			<h1>Create table:</h1>
			<table border="0" width="800px" cellspacing="0" cellpadding="0px">
				<caption>Enter the required info to make yourself a table in MS Excel format</caption>
				<tr>
					<td class="label">
						Choose your field delimitor:
					</td>
					<td class="field">
					<table border="0">
						<tr>
						<td>
						<tr><td><input type="radio" name="field_delimitor" value='comma' >comma(,)</td></tr> 
						<tr><td><input type="radio" name="field_delimitor" value='tab' checked>tab(well, a tab...)</td></tr>
						</tr>
					</table>
					</td>
					<td class="label3">
						Choose your decimal sign:
					</td>
					<td class="field2">
						<table border="0">
						<tr>
						<td>
						<tr><td><input type="radio" name="decimal" value='comma' >comma(,)</td></tr> 
						<tr><td><input type="radio" name="decimal" value='dot' checked>dot(.)<td></tr>
						</td>
						</tr>
						</table>
					</td>
					<td class="label3">
						For protein codinng genes</br>
						- include individual codon positions?
					</td>
					<td class="field2">
						<table border="0">
						<tr>
						<td>
						<tr><td><input type="radio" name="codpos" value='yes' >yes</td></tr>
						<tr><td><input type="radio" name="codpos" value='no' checked>no</td></tr>
						</tr>
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
						if (isset($taxonsets)) {
							foreach ($taxonsets as $taxonset){ echo "<option value=\"$taxonset\">$taxonset</option> ";}
						}
						?>
						</select>
					</td>
				</tr>
				<tr>
					<td class="label4">
						...and/or a <br />list of codes:,
						:
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
						&nbsp;
					</td>
					<td class="label4">
						Check to select</br> your alignment/gene:
					</td>
					<td class="field1">
						<?php $i = 0;
								echo "<table border=\"0\">";
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
					<td> &nbsp; </td>
					<td class="field1">
						<textarea name="codes" rows="14"></textarea>
					</td><td></td>
					<td colspan=3 style="font-size:16px">Keep in mind that this might take</br> a few minutes for large data sets!
				</br>
				</br>
					<input type="submit" name="make_gene_table" value="Create gene table" />
				</td>
			</table>
		</td>
	</tr>
	<tr>
		
	</tr>
</table>

</form>

</div> <!-- end content -->

<?php
make_footer($date_timezone, $config_sitename, $version, $base_url);
?>
</body>
</html>
