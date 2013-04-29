<?php
// #################################################################################
// #################################################################################
// Voseq admin/add_gene.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: add and update genes and their information:
//  - reading frames
//  - number of basepairs 
//  - description
//
// #################################################################################


// #################################################################################
// Section: include functions
// #################################################################################

//check admin login session
include'../login/auth-admin.php';

error_reporting (E_ALL ^ E_NOTICE);

// includes
#include '../login/redirect.html';
ob_start();//Hook output buffer - disallows web printing of file info...
include'../conf.php';
ob_end_clean();//Clear output buffer//includes
include 'adfunctions.php'; // administrator functions
include 'admarkup-functions.php';
include '../includes/validate_coords.php';

// process title
$title = $config_sitename;

// need dojo?
$dojo = true;

// which dojo?
$whichDojo[] = 'Tooltip';
$whichDojo[] = 'ComboBox';

// to indicate this is an administrator page
$admin = true;





// #################################################################################
// form not yet submitted
// display initial form and empty
// #################################################################################
if ($_GET['new'] || $_POST['submitNewIntrons'])
	// brand new record or in creation
	{
	
	// capture values from record in creation
	if ($_POST['submitNewIntrons']){
		$geneCode      = trim($_POST['geneCode']);
		$description = utf8_encode (trim($_POST['description']));
		$length     = trim($_POST['length']);
		$readingframe     = $_POST['readingframe'];
		$notes     = $_POST['notes'];
		$genetype = $_POST['genetype'];
		$genetype_other_name = trim($_POST['genetype_other_name']);
		$numIntrons = $_POST['numIntrons'];
		$protcoding = $_POST['protcoding'];
		$aligned = $_POST['aligned'];
		$genetic_code = $_POST['genetic_code'];
		if(isset($_POST['intron_start'])) { 
			$i = 1; 
			foreach ($_POST['intron_start'] as $is) {
			$intron_start[$i] = $is; $i++;
			}
		}
		if(isset($_POST['intron_end'])) { 
			$i = 1; foreach ($_POST['intron_end'] as $ie) {
			$intron_end[$i] = $ie; $i++;
			}
		}
	}
	else {
		$protcoding = 'no';
		$genetype = 'other';
		$numIntrons = '0';
		$aligned = 'no';
		$genetic_code = '0';
	}
	// print html headers
	include_once('../includes/header.php');

	// print navegation bar
	admin_nav();

	// begin HTML page content
	echo "<div id=\"content\">";
	?>
	

<table border="0" width="960px"> <!-- super table -->
<tr><td valign="top">
	<form action="<?php echo $_SERVER['PHP_SELF']; ?>" method="post">

<b>Create a definition for your alignment/gene. Specify "Reading frame" if you want to create datasets by codon positions.</b>

<table width="800" border="0"> <!-- big parent table -->
<tr><td valign="top">
	<table border="0" cellspacing="10"> <!-- table child 1 -->
	<tr><td>
	<table width="575" cellspacing="0" border="0">
	<caption>Gene information</caption>
		<tr>
			<td class="label">Gene code</td>
			<td class="field">
				<input size="12" maxlength="40" type="text" name="geneCode" value="<?php echo $geneCode; ?>"/>
				</select></td>
			<td class="label3">Aligned or not?</td>
			<td class="field2" colspan="1">
				<input type="radio" name="aligned" value="yes" <?php if ($aligned == "yes") { echo " checked "; }?> />yes 
				<input type="radio" name="aligned" value="no" <?php if ($aligned == "no") { echo " checked "; }?> >no
			</td>
			<td class="label3">Length (if aligned):</td>
			<td class="field2">
				<input size="10" maxlength="40" type="text" name="length" value="<?php echo $length; ?>"/>
				</select></td>
		</tr>
		<tr>
			<td class="label">Description</td>
			<td class="field" colspan = "5">
					<input size="80" maxlength="500" type="text" name="description" value="<?php echo $description; ?>"/>
				</select></td>
		</tr>
		<tr>
			<td class="label">Notes:
			</td>
			<td class="field" colspan="5">
					<input size="80" maxlength="500" type="text" name="notes" value="<?php echo $notes; ?>"/>
			</td>
		</tr>
		<tr>
			<td class="label">Gene type:</td>
			<td class="field" colspan = "5">
				<input type="radio" name="genetype" value="mitochondrial" <?php if ($genetype == "mitochondrial") { echo " checked "; }?> >mitochondrial 
				<input type="radio" name="genetype" value="nuclear" <?php if ($genetype == "nuclear") { echo " checked "; }?> >nuclear 
				<input type="radio" name="genetype" value="ribosomal" <?php if ($genetype == "ribosomal") { echo " checked "; }?> >ribosomal
				<input type="radio" name="genetype" value="plastid" <?php if ($genetype == "plastid") { echo " checked "; }?> >plastid
				<input type="radio" name="genetype" value="other" <?php if ($genetype == "other") { echo " checked "; }?> >other:
				<input size="16" maxlength="80" type="text" name="genetype_other_name" />
				</select>
			</td>
		</tr>
		</table>
		</td>
		</tr>
		<tr>
		<td>
		<table width="575" cellspacing="0" border="0">
		<tr><caption colspan=6>Protein coding genes</caption></tr>
		<tr><td class="label"></td><td class="field" colspan="6">This section will only be applied if 'Protein coding' is set to 'yes'</td></tr>
		<tr>
			<td class="label">Protein coding</td>
			<td class="field" colspan="2">
				<input type="radio" name="protcoding" value="yes" <?php if ($protcoding == "yes") { echo " checked "; }?> />yes 
				<input type="radio" name="protcoding" value="no" <?php if ($protcoding == "no") { echo " checked "; }?> >no
			</td>
			<td class="label3">Reading frame</td>
			<td class="field2" colspan="2">
				<input type="radio" name="readingframe" value="1" <?php if ($readingframe == "1") { echo " checked "; }?>>1 
				<input type="radio" name="readingframe" value="2" <?php if ($readingframe == "2") { echo " checked "; }?>>2 
				<input type="radio" name="readingframe" value="3" <?php if ($readingframe == "3") { echo " checked "; }?>>3
			</td>
		</tr>
		<tr>
			<td class="label">Translation</br>table</td>
			<td class="field" colspan = "5">
					<select name="genetic_code" size="1" style=" BORDER-BOTTOM: outset; BORDER-LEFT: 
			outset; BORDER-RIGHT: outset; BORDER-TOP: outset; FONT-FAMILY: 
			Arial; FONT-SIZE: 12px"> 
			  <!-- create a pulldown-list with all taxon set names in the db -->
			    <option value=0 <?php if ($genetic_code == "0") { echo "selected "; }?>>
				<option value=1 <?php if ($genetic_code == "1") { echo "selected "; }?>>Standard
                <option value=2 <?php if ($genetic_code == "2") { echo "selected "; }?>>Vertebrate Mitochondrial
                <option value=3 <?php if ($genetic_code == "3") { echo "selected "; }?>>Yeast Mitochondrial
                <option value=4 <?php if ($genetic_code == "4") { echo "selected "; }?>>Mold, Protozoan and Coelenterate Mitochondrial. Mycoplasma, Spiroplasma
                <option value=5 <?php if ($genetic_code == "5") { echo "selected "; }?>>Invertebrate Mitochondrial
                <option value=6 <?php if ($genetic_code == "6") { echo "selected "; }?>>Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
                <option value=9 <?php if ($genetic_code == "9") { echo "selected "; }?>>Echinoderm Mitochondrial
                <option value=10 <?php if ($genetic_code == "10") { echo "selected "; }?>>Euplotid Nuclear
                <option value=11 <?php if ($genetic_code == "11") { echo "selected "; }?>>Bacterial and Plant Plastid
                <option value=12 <?php if ($genetic_code == "12") { echo "selected "; }?>>Alternative Yeast Nuclear
                <option value=13 <?php if ($genetic_code == "13") { echo "selected "; }?>>Ascidian Mitochondrial
                <option value=14 <?php if ($genetic_code == "14") { echo "selected "; }?>>Flatworm Mitochondrial
                <option value=15 <?php if ($genetic_code == "15") { echo "selected "; }?>>Blepharisma Macronuclear
                <option value=16 <?php if ($genetic_code == "16") { echo "selected "; }?>>Chlorophycean Mitochondrial
                <option value=21 <?php if ($genetic_code == "21") { echo "selected "; }?>>Trematode Mitochondrial
                <option value=22 <?php if ($genetic_code == "22") { echo "selected "; }?>>Scenedesmus obliquus mitochondrial
                <option value=23 <?php if ($genetic_code == "23") { echo "selected "; }?>>Thraustochytrium mitochondrial code
				</select>
				</select></td>
		</tr>
		</table>
		</td>
		</tr>
		<tr>
		<td>
		<table width="575" cellspacing="0" border="0">
		<tr><caption colspan=6>Introns</caption></tr>
		<tr>
			<td class="label">Introns</td>
			<td class="field" colspan="5">
				Do you have introns in your alignment?</br> Enter number and press "Update/Add introns"</br>
				Then enter starting and ending nucelotide number for each intron region.</br>
				To remove introns just select fewer and press button again.</br>
				<!-- <form id="formName" action="<?php //echo $_SERVER['PHP_SELF'];?>" method="get"> -->
				<input size="3" maxlength="10" type="number" name="numIntrons" min="0" max="100" value="<?php echo $numIntrons; ?>"/>
				<input type="submit" name="submitNewIntrons" value="Update/Add introns" />
				<!-- </form> -->
			</td>
		</tr>
		<?php
		if($_POST['submitNewIntrons'] && $numIntrons > 0 ){
			$outp =  "";
			for ($i = 1; $i <= $numIntrons; $i++) {
				if ($i % 2 != 0) { $outp .= "</tr><td class=\"label\">";}
				else { $outp .= "<td class=\"label3\">";}
				$outp .= "Intron $i:</td>";
				if ($i % 2 != 0) { $outp .= "<td class=\"field\"";}
				else { $outp .= "<td class=\"field2\"";}
				$outp .=    " colspan=\"2\">
							<input size=\"4\" maxlength=\"5\" type=\"text\" name=\"intron_start[$i]\" value=\"$intron_start[$i]\"/>
							-
							<input size=\"4\" maxlength=\"5\" type=\"text\" name=\"intron_end[$i]\" value=\"$intron_end[$i]\"/>
							</td>";
				if ($i % 2 == 0) { $outp .= "</tr>";}
			}
			echo $outp;
		}
		?>
		</tr>
		<tr>
			<td></td><td></td><td></td><td></td><td></td><td>
				<input type="submit" name="submitNew" value="Add gene" />
			</td>
		</tr>
	</table>
	
	</td></tr>
	</table><!-- end table child 2 -->

</td></tr>
</table><!-- end big parent table -->

</td>
<td class="sidebar" valign="top">
	<?php admin_make_sidebar(); ?>
</td>
</tr>
</table> <!-- end super table -->

</form>
</div> <!-- end content -->

<!-- standard page footer begins -->
<?php
make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
?>
	
	<?php
	}
elseif ($_POST['submitNew']) {
	// set up error list array
	$errorList = array();
	
	//validate text input fields
	if (trim($_POST['geneCode']) == '')
		{
		$errorList[] = "Invalid entry: <b>Gene code</b></br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You must specify a gene code to proceed!";
		}

	$geneCode      = trim($_POST['geneCode']);
	$description = utf8_encode (trim($_POST['description']));
	$length     = trim($_POST['length']);
	$readingframe     = $_POST['readingframe'];
	$notes     = $_POST['notes'];
	$genetype = $_POST['genetype'];
	$genetype_other_name = trim($_POST['genetype_other_name']);
	$numIntrons = $_POST['numIntrons'];
	$protcoding = $_POST['protcoding'];
	$aligned = $_POST['aligned'];
	$genetic_code = $_POST['genetic_code'];
	if(isset($_POST['intron_start'])) { 
		$i = 1; 
		foreach ($_POST['intron_start'] as $is) {
		$intron_start[$i] = $is; $i++;
		}
	}
	if(isset($_POST['intron_end'])) { 
		$i = 1; foreach ($_POST['intron_end'] as $ie) {
		$intron_end[$i] = $ie; $i++;
		}
	}
	// fixing genetype and other
	if ($genetype == "other" && isset($genetype_other_name)){
		if ($genetype_other_name != '') {	
			$genetype == $genetype_other_name;
		}
	}
	
	// check if values are ok
	if (isset($length) && $aligned == 'yes'){
		if (is_numeric($length) && !$length == ''){$ltest = $length; $length = "'$length'";  }
		else {
			$errorList[] = "Invalid entry: <b>Length = \"$length\"</b></br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Length needs to be an integer!";
			}
	}
	else {$length = 'NULL'; $ltest = 'NULL';}

	// delete values if not prot coding or check values if it is
	if ($protcoding == "no") {
		$genetic_code = "NULL";
		$readingframe = "NULL";
	}
	else {
		if ($genetype == 'ribosomal'){
			$errorList[] = "Invalid entry: <b>Ribosomal genetype AND protein coding has been chosen</b></br>
					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It doesnt make sense!!";
		}
		if ($readingframe != '1' && $readingframe != '2' && $readingframe != '3'){
			$errorList[] = "Invalid entry: <b>Protein coding has been chosen but no reading frame is set!</b></br>
					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
		}
		$genetic_codes_arr = array('1','2','3','4','5','6','9','10','11','12','13','14','15','16','21','22','23');
		$genetic_codes_mito = array('2','3','4','5','9','13','14','16','21','22','23');
		$genetic_codes_nuc = array('1','6','10','12','15');
		if (!in_array($genetic_code, $genetic_codes_arr)){
			$errorList[] = "Invalid entry: <b>Protein coding has been chosen but no genetic code is set!</b></br>
					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
		}
		else {
			if ($genetype == 'mitochondrial' && !in_array($genetic_code, $genetic_codes_mito)){
				$errorList[] = "Invalid entry: <b>Mitochondrial genetype has been chosen but no mitochondrial genetic code!</b></br>
						&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
			}
			if ($genetype == 'nuclear' && !in_array($genetic_code, $genetic_codes_nuc)){
				$errorList[] = "Invalid entry: <b>Nuclear gene type has been chosen but no nuclear genetic code!</b></br>
						&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
			}
			if ($genetype == 'plastid' && $genetic_code != '11' ){
				$errorList[] = "Invalid entry: <b>Plastid gene type has been chosen but no plastid genetic code!</b></br>
						&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
			}
		}
	}

	// test introns
	if ($length == 'NULL' && $numIntrons > 0) {
		$errorList[] = "Invalid entry: <b>Introns have been speciefied but alignment length</br> (or 'aligned=yes') is not set</b></br>
				&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Go back and edit!!";
	}
	if ($numIntrons > 0 && !isset($intron_start)){
					$errorList[] = "Invalid entry: <b>$numIntrons introns are speciefied</br>but starting nucleotide number(s) is not set</b></br>
							&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Set introns to 0 or specify starting point(s)!";
	}
	if ($numIntrons > 0 && !isset($intron_start)){
					$errorList[] = "Invalid entry: <b>$numIntrons introns are speciefied</br>but ending nucleotide number(s) is not set</b></br>
							&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Set introns to 0 or specify ending point(s)!";
	}
	if (isset($intron_start) && isset($intron_end)){
		for ($i = 1; $i <= $numIntrons; $i++) {
			if (!is_numeric($intron_start[$i]) || !is_numeric($intron_end[$i])){
				if (!is_numeric($intron_start[$i])){
						$errorList[] = "Invalid entry: <b>Intron $i start value = \"$intron_start[$i]\"</b></br>
										&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It needs to be an integer!";
					}
				if (!is_numeric($intron_end[$i])){
						$errorList[] = "Invalid entry: <b>Intron $i end value = \"$intron_end[$i]\"</b></br>
										&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It needs to be an integer!";
				}
			}
			else {
				if ($intron_start[$i] > $intron_end[$i]){ 
						$errorList[] = "Invalid entry: <b>Intron $i start value($intron_start[$i])</b> </br>
										is higher than its end value(\"$intron_end[$i]\")</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;PLease correct this!!";
				}
				if ($intron_start[$i] > $ltest || $intron_end[$i] > $ltest){
						$errorList[] = "Invalid entry: <b>Intron $i start value($intron_start[$i]) or end value($intron_end[$i])</b></br>
										is higher than total specified gene length</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It doesnt make sense!!";
				}
				if ($numIntrons > 1 && $i > 1){
					if ($intron_start[$i] <= $intron_end[$i-1]){
							$errorList[] = "Invalid entry: <b>Intron $i start value($intron_start[$i]) is the same or lower than Intron ". ($i-1) ."(".$intron_end[$i-1].")</b></br>
											</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Introns should be speciefied in correct order!!";
					}
				}
			}
		}
	}
	// check for errors
	// if none found ...
	if (sizeof($errorList) == 0 )
		{
		// open database connection
		$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
		//select database
		mysql_select_db($db) or die ('Unable to content');
		if( function_exists(mysql_set_charset) ) {
			mysql_set_charset("utf8");
		}
		
		// check for duplicate geneCode
		$querygCode = "SELECT * FROM ". $p_ . "genes WHERE geneCode='$geneCode'";
		$resultgCode = mysql_query($querygCode) or die ("Error in query: $querygCode. " . mysql_error());
		if (mysql_num_rows($resultgCode) > 0)
			{
			// process title
			$title = "$config_sitename - Error, duplicate gene code";
			
			// print html headers
			include_once('../includes/header.php');
			admin_nav();
			
			// begin HTML page content
			echo "<div id=\"content_narrow\">";
			echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
					<tr><td valign=\"top\">";
			echo "<img src=\"../images/warning.png\" alt=\"\">
						The record's <b>code</b> you entered is already preoccupied.<br />There can't be two genes with the same gene code!.<br />Please click \"Go back\" in your browser and enter a different gene code.</br></br></br>";
		echo "<form action=\"" . $_SERVER['PHP_SELF'] . "\" method=\"post\">";
		foreach ($_POST as  $key=>$value) { 
			if ($key == "intron_start" || $key == "intron_end") {
					for ($i = 1; $i <= $numIntrons; $i++){
						echo "<input type=\"hidden\" name=\"intron_start[]\" value=\"$intron_start[$i]\">" ;
						echo "<input type=\"hidden\" name=\"intron_end[]\" value=\"$intron_end[$i]\">" ;
					}
			}
			else {
				echo "<input type=\"hidden\" name=\"$key\" value=\"$value\">" ;
			}
		}
		echo "<input type=\"submit\" name=\"submitNewIntrons\" value=\"Go Back And Edit Here\" /></td></form>";
			echo "<td class=\"sidebar\" valign=\"top\">";
				admin_make_sidebar();
			echo "</td>";
			echo "</tr>
					</table> <!-- end super table -->
					</div> <!-- end content -->";
			make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
			echo "\n</body>\n</html>";
			exit();
			}
		else
			{	
			//setting the edits add values
			// $editsadd = "Added by ". $_SESSION['SESS_FIRST_NAME']. " ". $_SESSION['SESS_LAST_NAME'] ." on ";
			// mysql_query("time for add-list");
			// $querytime = "SELECT NOW()";
			// $resulttime = mysql_query($querytime) or die ("Error in query: $querytime. " . mysql_error());
			// $rowtime    = mysql_result($resulttime,0);
			// $editsadd = $editsadd . $rowtime;
			// mysql_query("set names utf8");
			
			// fixing intron values into string
			if ($numIntrons > 0){
				for ($i = 1; $i <= $numIntrons; $i++){
					$introns[] = "$intron_start[$i]-$intron_end[$i]";
				}
				$intron_string = implode(";", $introns);
				$intron_string = "'$intron_string'";
			}
			else { $intron_string = "NULL"; }

			
			// generate and execute query
			$gquery_insert = "INSERT INTO ". $p_ . "genes(geneCode, length, description, readingframe, notes, genetype, prot_code, intron, genetic_code, aligned, timestamp)"; 
			//$gquery_values = " VALUES ('$geneCode', '$length', '$description', '$readingframe', '$notes','$genetype', '$protcoding', ";
			$gquery_values = " VALUES ('$geneCode', $length, '$description', $readingframe,'$notes','$genetype','$protcoding', $intron_string, $genetic_code,'$aligned', NOW())";
	
			// if 	($numIntrons > 0){	$gquery .= "'$intron_string', ";}
			// else {$gquery .= "NULL, ";}
			// if 	($genetic_code > 0){	$gquery .= "'$genetic_code', ";}
			// else {$gquery .= "NULL, ";}
			// $gquery .= "'$aligned', NOW())";
			$gquery = "$gquery_insert $gquery_values";
			$gresult = mysql_query($gquery) or die ("Error in query: $query. " . mysql_error());
			
			// process title
			$title = "$config_sitename - Gene " . $geneCode . " created";

			// print html headers
			include_once('../includes/header.php');

			// print navegation bar
			admin_nav();

			// begin HTML page content
			echo "<div id=\"content_narrow\">";
			echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
					<tr><td valign=\"top\">";
			// print result
			echo "<span class=\"title\"><img src=\"images/success.png\" alt=\"\"> Gene creation was successful!</span>";
			}

			echo "<td class=\"sidebar\" valign=\"top\">";
				admin_make_sidebar();
			echo "</td>";
			echo "</tr>
				</table> <!-- end super table -->
				</div> <!-- end content -->";
			make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
		}
	else
		{
		// error found
		
		// get title
		$title = "$config_sitename - Error, missing info";
		
		// print html headers
		include_once('../includes/header.php');
		admin_nav();
		
		// begin HTML page content
		echo "<div id=\"content_narrow\">";
		echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
				<tr><td valign=\"top\">";
		echo "<img src=\"../images/warning.png\" alt=\"\"> The following errors were encountered:";
		echo '<br>';
		echo '<ul>';
		for ($x=0; $x<sizeof($errorList); $x++)
			{
			echo "<li>$errorList[$x]";
			}
		echo "<form action=\"" . $_SERVER['PHP_SELF'] . "\" method=\"post\">";
		foreach ($_POST as  $key=>$value) { 
			if ($key == "intron_start" || $key == "intron_end") {
					for ($i = 1; $i <= $numIntrons; $i++){
						echo "<input type=\"hidden\" name=\"intron_start[]\" value=\"$intron_start[$i]\">" ;
						echo "<input type=\"hidden\" name=\"intron_end[]\" value=\"$intron_end[$i]\">" ;
					}
			}
			else {
				echo "<input type=\"hidden\" name=\"$key\" value=\"$value\">" ;
			}
		}
		echo "</br><input type=\"submit\" name=\"submitNewIntrons\" value=\"Go Back And Edit Here\" /></form>";
		echo "</ul></td>";
		echo "<td class=\"sidebar\" valign=\"top\">";
			admin_make_sidebar();
		echo "</td>";
		echo "</tr>
				</table> <!-- end super table -->
				</div> <!-- end content -->";
		make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
		}
	}

// #################################################################################
// record to update
// get values to prefill fields
// #################################################################################
elseif (!$_POST['submitNewIntrons'] && !$_POST['submitNoNew'] && $_GET['geneCode'] || $_POST['updateNewIntrons']) {
	$geneCode1 = $_GET['geneCode'];
	@$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
	//select database
	mysql_select_db($db) or die ('Unable to content');
	if( function_exists(mysql_set_charset) ) {
		mysql_set_charset("utf8");
	}

	// check for duplicate code
	$query1  = "SELECT id, geneCode, length, description, readingframe, notes, intron, genetype, prot_code, aligned, genetic_code
				FROM ". $p_ . "genes WHERE geneCode='$geneCode1'";
	$result1 = mysql_query($query1) or die ("Error in query: $query1. " . mysql_error());
	$row1    = mysql_fetch_object($result1);
	
	// if update/add introns used store and output variables again
	if ($_POST['updateNewIntrons']){
		$geneCode      = trim($_POST['geneCode']);
		$id = $_POST['id'];
		$description = utf8_encode (trim($_POST['description']));
		$length     = trim($_POST['length']);
		$readingframe     = $_POST['readingframe'];
		$notes     = $_POST['notes'];
		$genetype = $_POST['genetype'];
		$genetype_other_name = trim($_POST['genetype_other_name']);
		$numIntrons = $_POST['numIntrons'];
		$protcoding = $_POST['protcoding'];
		$aligned = $_POST['aligned'];
		$genetic_code = $_POST['genetic_code'];
		if(isset($_POST['intron_start'])) { 
			$i = 1; 
			foreach ($_POST['intron_start'] as $is) {
			$intron_start[$i] = $is; $i++;
			}
		}
		if(isset($_POST['intron_end'])) { 
			$i = 1; foreach ($_POST['intron_end'] as $ie) {
			$intron_end[$i] = $ie; $i++;
			}
		}
	}
	else {
		$geneCode      = $row1->geneCode;
		$description = utf8_decode($row1->description);
		$length     = $row1->length;
		$readingframe     = $row1->readingframe;
		$notes     = $row1->notes;
		$protcoding = $row1->prot_code;
		if ($protcoding == 'notset'){$protcoding = 'no';}
		$genetype = $row1->genetype;
		$aligned = $row1->aligned;
		if ($aligned == 'notset'){$aligned = 'no';}
		$genetic_code = $row1->genetic_code;
		$id = $row1->id;
		if ($row1->intron != ''){
			$introns_raw = explode(";", $row1->intron);
					$numIntrons = count($introns_raw);
			for ($i = 1; $i <= $numIntrons; $i++){
				$j = $i-1;
				$intron_raw = explode("-", $introns_raw[$j]);
				$intron_start[$i] = $intron_raw[0];
				$intron_end[$i] = $intron_raw[1];
			}
		}
		else { $numIntrons = 0; }
		
		
	}
	// get title
	$title = "$config_sitename - Edit " . $geneCode1;
				
	// print html headers
	include_once('../includes/header.php');
	admin_nav();
				
	// begin HTML page content
	echo "<div id=\"content\">";
	
	// Delete button
	//echo "<button class='delete' style='background-color:red;color:white' id='delete_gene' name='" . $geneCode . "' class='delete'>Delete me</button>";
?>
<table border="0" width="960px"> <!-- super table -->
<tr><td valign="top">	


<table width="800" border="0"> <!-- big parent table -->
<?php 	// Delete button
echo "<tr><td><button id='delete_gene' name='" . $geneCode;
echo "'>Delete me</button>";
echo "</td></tr>";
?>
<form action="<?php echo $_SERVER['PHP_SELF']; ?>" method="post">
<tr><td valign="top">
	<table border="0" cellspacing="10"> <!-- table child 1 -->
	<tr><td>
		<!-- 	input id of this record also, useful for changing the code -->
		<input type="hidden" name="id" value="<?php echo $id; ?>" />
		<!-- 	end input id -->
	<table width="575" cellspacing="0" border="0">
	<caption>Gene information</caption>
		<tr>
			<td class="label">Gene code</td>
			<td class="field">
				<input size="12" maxlength="40" type="text" name="geneCode" value="<?php echo $geneCode; ?>"/>
				</select></td>
			<td class="label3">Aligned or not?</td>
			<td class="field2" colspan="1">
				<input type="radio" name="aligned" value="yes" <?php if ($aligned == "yes") { echo " checked "; }?> />yes 
				<input type="radio" name="aligned" value="no" <?php if ($aligned == "no") { echo " checked "; }?> >no
			</td>
			<td class="label3">Length (if aligned):</td>
			<td class="field2">
				<input size="10" maxlength="40" type="text" name="length" value="<?php echo $length; ?>"/>
				</select></td>
		</tr>
		<tr>
			<td class="label">Description</td>
			<td class="field" colspan = "5">
					<input size="80" maxlength="500" type="text" name="description" value="<?php echo $description; ?>"/>
				</select></td>
		</tr>
		<tr>
			<td class="label">Notes:
			</td>
			<td class="field" colspan="5">
					<input size="80" maxlength="500" type="text" name="notes" value="<?php echo $notes; ?>"/>
			</td>
		</tr>
		<tr>
			<td class="label">Gene type:</td>
			<td class="field" colspan = "5">
				<input type="radio" name="genetype" value="mitochondrial" <?php if ($genetype == "mitochondrial") { echo " checked "; }?> >mitochondrial 
				<input type="radio" name="genetype" value="nuclear" <?php if ($genetype == "nuclear") { echo " checked "; }?> >nuclear 
				<input type="radio" name="genetype" value="ribosomal" <?php if ($genetype == "ribosomal") { echo " checked "; }?> >ribosomal
				<input type="radio" name="genetype" value="plastid" <?php if ($genetype == "plastid") { echo " checked "; }?> >plastid
				<input type="radio" name="genetype" value="other" <?php if ($genetype == "other") { echo " checked "; }?> >other:
				<input size="16" maxlength="80" type="text" name="genetype_other_name" />
				</select>
			</td>
		</tr>
		</table>
		</td>
		</tr>
		<tr>
		<td>
		<table width="575" cellspacing="0" border="0">
		<tr><caption colspan=6>Protein coding genes</caption></tr>
		<tr><td class="label"></td><td class="field" colspan="6">This section will only be applied if 'Protein coding' is set to 'yes'</td></tr>
		<tr>
			<td class="label">Protein coding</td>
			<td class="field" colspan="2">
				<input type="radio" name="protcoding" value="yes" <?php if ($protcoding == "yes") { echo " checked "; }?> />yes 
				<input type="radio" name="protcoding" value="no" <?php if ($protcoding == "no") { echo " checked "; }?> >no
			</td>
			<td class="label3">Reading frame</td>
			<td class="field2" colspan="2">
				<input type="radio" name="readingframe" value="1" <?php if ($readingframe == "1") { echo " checked "; }?>>1 
				<input type="radio" name="readingframe" value="2" <?php if ($readingframe == "2") { echo " checked "; }?>>2 
				<input type="radio" name="readingframe" value="3" <?php if ($readingframe == "3") { echo " checked "; }?>>3
			</td>
		</tr>
		<tr>
			<td class="label">Translation</br>table</td>
			<td class="field" colspan = "5">
					<select name="genetic_code" size="1" style=" BORDER-BOTTOM: outset; BORDER-LEFT: 
			outset; BORDER-RIGHT: outset; BORDER-TOP: outset; FONT-FAMILY: 
			Arial; FONT-SIZE: 12px"> 
			  <!-- create a pulldown-list with all taxon set names in the db -->
			    <option value=0 <?php if ($genetic_code == "0") { echo "selected "; }?>>
				<option value=1 <?php if ($genetic_code == "1") { echo "selected "; }?>>Standard
                <option value=2 <?php if ($genetic_code == "2") { echo "selected "; }?>>Vertebrate Mitochondrial
                <option value=3 <?php if ($genetic_code == "3") { echo "selected "; }?>>Yeast Mitochondrial
                <option value=4 <?php if ($genetic_code == "4") { echo "selected "; }?>>Mold, Protozoan and Coelenterate Mitochondrial. Mycoplasma, Spiroplasma
                <option value=5 <?php if ($genetic_code == "5") { echo "selected "; }?>>Invertebrate Mitochondrial
                <option value=6 <?php if ($genetic_code == "6") { echo "selected "; }?>>Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
                <option value=9 <?php if ($genetic_code == "9") { echo "selected "; }?>>Echinoderm Mitochondrial
                <option value=10 <?php if ($genetic_code == "10") { echo "selected "; }?>>Euplotid Nuclear
                <option value=11 <?php if ($genetic_code == "11") { echo "selected "; }?>>Bacterial and Plant Plastid
                <option value=12 <?php if ($genetic_code == "12") { echo "selected "; }?>>Alternative Yeast Nuclear
                <option value=13 <?php if ($genetic_code == "13") { echo "selected "; }?>>Ascidian Mitochondrial
                <option value=14 <?php if ($genetic_code == "14") { echo "selected "; }?>>Flatworm Mitochondrial
                <option value=15 <?php if ($genetic_code == "15") { echo "selected "; }?>>Blepharisma Macronuclear
                <option value=16 <?php if ($genetic_code == "16") { echo "selected "; }?>>Chlorophycean Mitochondrial
                <option value=21 <?php if ($genetic_code == "21") { echo "selected "; }?>>Trematode Mitochondrial
                <option value=22 <?php if ($genetic_code == "22") { echo "selected "; }?>>Scenedesmus obliquus mitochondrial
                <option value=23 <?php if ($genetic_code == "23") { echo "selected "; }?>>Thraustochytrium mitochondrial code
				</select>
				</select></td>
		</tr>
		</table>
		</td>
		</tr>
		<tr>
		<td>
		<table width="575" cellspacing="0" border="0">
		<tr><caption colspan=6>Introns</caption></tr>
		<tr>
			<td class="label">Introns</td>
			<td class="field" colspan="5">
				Do you have introns in your alignment?</br> Enter number and press "Update/Add introns"</br>
				Then enter starting and ending nucelotide number for each intron region.</br>
				To remove introns just select fewer and press button again.</br>
				<!-- <form id="formName" action="<?php //echo $_SERVER['PHP_SELF'];?>" method="get"> -->
				<input size="3" maxlength="10" type="number" name="numIntrons" min="0" max="100" value="<?php echo $numIntrons; ?>"/>
				<input type="submit" name="updateNewIntrons" value="Update/Add introns" />
				<!-- </form> -->
			</td>
		</tr>
		<?php
		if($_POST['updateNewIntrons'] || $numIntrons > 0 ){
			$outp =  "";
			for ($i = 1; $i <= $numIntrons; $i++) {
				if ($i % 2 != 0) { $outp .= "</tr><td class=\"label\">";}
				else { $outp .= "<td class=\"label3\">";}
				$outp .= "Intron $i:</td>";
				if ($i % 2 != 0) { $outp .= "<td class=\"field\"";}
				else { $outp .= "<td class=\"field2\"";}
				$outp .=    " colspan=\"2\">
							<input size=\"4\" maxlength=\"5\" type=\"text\" name=\"intron_start[$i]\" value=\"$intron_start[$i]\"/>
							-
							<input size=\"4\" maxlength=\"5\" type=\"text\" name=\"intron_end[$i]\" value=\"$intron_end[$i]\"/>
							</td>";
				if ($i % 2 == 0) { $outp .= "</tr>";}
			}
			echo $outp;
		}
		?>
		</tr>
		<tr>
					<td></td><td></td><td></td><td></td><td></td><td>
				<input type="submit" name="submitNoNew" value="Update gene" />
			</td>
		</tr>
	</table>
	
	</td></tr>
	</table><!-- end table child 2 -->

</td></tr>
</table><!-- end big parent table -->

</td>

<td class="sidebar" valign="top">
	<?php admin_make_sidebar();  ?>
</td>

</tr>
</table>
</table> <!-- end super table -->

</form>
</div> <!-- end content -->

	<?php
	// close database connection
	mysql_close($connection);

	make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);

}
elseif ($_POST['submitNoNew']) {
	// set up error list array
	$errorList = array();
	
	//validate text input fields
	if (trim($_POST['geneCode']) == '')
		{
		$errorList[] = "Invalid entry: <b>Gene code</b></br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You must specify a gene code to proceed!";
		}
	$id1       = $_POST['id'];
	$geneCode1      = trim($_POST['geneCode']);
	$description = utf8_encode (trim($_POST['description']));
	$length     = trim($_POST['length']);
	$readingframe     = $_POST['readingframe'];
	$notes     = $_POST['notes'];
	$genetype = $_POST['genetype'];
	$genetype_other_name = trim($_POST['genetype_other_name']);
	$numIntrons = $_POST['numIntrons'];
	$protcoding = $_POST['protcoding'];
	$aligned = $_POST['aligned'];
	$genetic_code = $_POST['genetic_code'];
	if(isset($_POST['intron_start'])) { 
		$i = 1; 
		foreach ($_POST['intron_start'] as $is) {
		$intron_start[$i] = $is; $i++;
		}
	}
	if(isset($_POST['intron_end'])) { 
		$i = 1; foreach ($_POST['intron_end'] as $ie) {
		$intron_end[$i] = $ie; $i++;
		}
	}
	// fixing genetype and other
	if ($genetype == "other" && isset($genetype_other_name)){
		if ($genetype_other_name != '') {	
			$genetype == $genetype_other_name;
		}
	}
	
	// check if values are ok
	if (isset($length) && $aligned == 'yes'){
		if (is_numeric($length) && !$length == ''){$ltest = $length; $length = "'$length'";  }
		else {
			$errorList[] = "Invalid entry: <b>Length = \"$length\"</b></br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Length needs to be an integer!";
			}
	}
	else {$length = 'NULL'; $ltest = 'NULL';}

	// delete values if not prot coding or check values if it is
	if ($protcoding == "no") {
		$genetic_code = "NULL";
		$readingframe = "NULL";
	}
	else {
		if ($genetype == 'ribosomal'){
			$errorList[] = "Invalid entry: <b>Ribosomal genetype AND protein coding has been chosen</b></br>
					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It doesnt make sense!!";
		}
		if ($readingframe != '1' && $readingframe != '2' && $readingframe != '3'){
			$errorList[] = "Invalid entry: <b>Protein coding has been chosen but no reading frame is set!</b></br>
					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
		}
		else {$readingframe = "'$readingframe'";}
		$genetic_codes_arr = array('1','2','3','4','5','6','9','10','11','12','13','14','15','16','21','22','23');
		$genetic_codes_mito = array('2','3','4','5','9','13','14','16','21','22','23');
		$genetic_codes_nuc = array('1','6','10','12','15');
		if (!in_array($genetic_code, $genetic_codes_arr)){
			$errorList[] = "Invalid entry: <b>Protein coding has been chosen but no genetic code is set!</b></br>
					&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
		}
		else {
			if ($genetype == 'mitochondrial' && !in_array($genetic_code, $genetic_codes_mito)){
				$errorList[] = "Invalid entry: <b>Mitochondrial genetype has been chosen but no mitochondrial genetic code!</b></br>
						&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
			}
			if ($genetype == 'nuclear' && !in_array($genetic_code, $genetic_codes_nuc)){
				$errorList[] = "Invalid entry: <b>Nuclear gene type has been chosen but no nuclear genetic code!</b></br>
						&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
			}
			if ($genetype == 'plastid' && $genetic_code != '11' ){
				$errorList[] = "Invalid entry: <b>Plastid gene type has been chosen but no plastid genetic code!</b></br>
						&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please go back and edit";
			}
		}
	}

	// test introns
	if ($length == 'NULL' && $numIntrons > 0) {
		$errorList[] = "Invalid entry: <b>Introns have been speciefied but alignment length</br> (or 'aligned=yes') is not set</b></br>
				&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Go back and edit!!";
	}
	if ($numIntrons > 0 && !isset($intron_start)){
					$errorList[] = "Invalid entry: <b>$numIntrons introns are speciefied</br>but starting nucleotide number(s) is not set</b></br>
							&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Set introns to 0 or specify starting point(s)!";
	}
	if ($numIntrons > 0 && !isset($intron_start)){
					$errorList[] = "Invalid entry: <b>$numIntrons introns are speciefied</br>but ending nucleotide number(s) is not set</b></br>
							&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Set introns to 0 or specify ending point(s)!";
	}
	if (isset($intron_start) && isset($intron_end)){
		for ($i = 1; $i <= $numIntrons; $i++) {
			if (!is_numeric($intron_start[$i]) || !is_numeric($intron_end[$i])){
				if (!is_numeric($intron_start[$i])){
						$errorList[] = "Invalid entry: <b>Intron $i start value = \"$intron_start[$i]\"</b></br>
										&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It needs to be an integer!";
					}
				if (!is_numeric($intron_end[$i])){
						$errorList[] = "Invalid entry: <b>Intron $i end value = \"$intron_end[$i]\"</b></br>
										&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It needs to be an integer!";
				}
			}
			else {
				if ($intron_start[$i] > $intron_end[$i]){ 
						$errorList[] = "Invalid entry: <b>Intron $i start value($intron_start[$i])</b> </br>
										is higher than its end value(\"$intron_end[$i]\")</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;PLease correct this!!";
				}
				if ($intron_start[$i] > $ltest || $intron_end[$i] > $ltest){
						$errorList[] = "Invalid entry: <b>Intron $i start value($intron_start[$i]) or end value($intron_end[$i])</br>
										is higher than total specified gene length ($ltest)</b></br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;It doesnt make sense!!";
				}
				if ($numIntrons > 1 && $i > 1){
					if ($intron_start[$i] <= $intron_end[$i-1]){
							$errorList[] = "Invalid entry: <b>Intron $i start value($intron_start[$i]) is the same or lower than Intron ". ($i-1) ."(".$intron_end[$i-1].")</b></br>
											</br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Introns should be speciefied in correct order!!";
					}
				}
			}
		}
	}

	// check for errors
	// if none found ...
	if (sizeof($errorList) == 0 )
		{
		// open database connection
		$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
		
		//select database
		mysql_select_db($db) or die ('Unable to content');
		if( function_exists(mysql_set_charset) ) {
			mysql_set_charset("utf8");
		}
		
		// check if submitted code is meant to replace old one
		// get old code
		$queryOldCode = "SELECT geneCode FROM ". $p_ . "genes WHERE id='$id1'";
		$resultOldCode = mysql_query($queryOldCode) or die ("Error in query: $queryOldCode. " . mysql_error());
		$rowOldCode    = mysql_fetch_object($resultOldCode);
		$oldCode = $rowOldCode->geneCode;
		// get new code
		$newCode = $geneCode1;
		//  if new code != old code
		if ($oldCode != $newCode) {
			// check for duplicate
			$queryCode1 = "SELECT geneCode FROM ". $p_ . "genes WHERE geneCode='$newCode'";
			$resultCode1 = mysql_query($queryCode1) or die ("Error in query: $queryCode1. " . mysql_error());
			if (mysql_num_rows($resultCode1) > 0)
				{
				// get title
				$title = "$config_sitename - Error, duplicate gene code";
				
				// print html headers
				include_once('../includes/header.php');
				admin_nav();
				
				// begin HTML page content
				echo "<div id=\"content_narrow\">";
				echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
						<tr><td valign=\"top\">";
				echo "<img src=\"../images/warning.png\" alt=\"\">
						The record's <b>gene code</b> ($newCode) you entered is already preoccupied.<br />
						There can't be two genes with the same gene code!.<br /><br />
						Please click \"Go back\" in your browser and enter a different code.</span>
						</td>";
				echo "<td class=\"sidebar\" valign=\"top\">";
				admin_make_sidebar(); 
				echo "</td>";
				echo "</tr>
					  </table> <!-- end super table -->
					  </div> <!-- end content -->";
				make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
				exit();
				}
			}
		// fixing intron values into string
			if ($numIntrons > 0){
				for ($i = 1; $i <= $numIntrons; $i++){
					$introns[] = "$intron_start[$i]-$intron_end[$i]";
				}
				$intron_string = implode(";", $introns);
				$intron_string = "'$intron_string'";
			}
			else { $intron_string = "NULL"; }
		// utf8 encode some fields
		$geneCode1 = $geneCode1;
		$description = utf8_encode($description);
		// generate and execute query UPDATE
		$query = "UPDATE ". $p_ . "genes SET geneCode='$geneCode1', length=$length, description='$description', readingframe=$readingframe, notes='$notes', 
				genetype='$genetype', prot_code='$protcoding',aligned='$aligned', intron=$intron_string, genetic_code=$genetic_code, timestamp=NOW() WHERE id='$id1'";
		$result = mysql_query($query) or die ("Error in query: $query. " . mysql_error());
		
		//update all sequences and primers with old genecode to new genecode
		if ($oldCode != $newCode){
			$tablelist = array("sequences", "primers");
			foreach ($tablelist as $tabLe){
				$querygC = "UPDATE ". $p_ . "$tabLe SET geneCode='$geneCode1' WHERE geneCode='$oldCode'";
				$resultgC = mysql_query($querygC) or die ("Error in query: $querygC. " . mysql_error());
			}
		}
		// get title
		$title = "$config_sitename - Record " . $geneCode1 . " updated";
				
		// print html headers
		include_once('../includes/header.php');
		admin_nav();
				
		// begin HTML page content
		echo "<div id=\"content_narrow\">";
		echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
				<tr><td valign=\"top\">";
		echo "<img src=\"images/success.png\" alt=\"\"> Record update was successful!";
		echo "</td>";
		echo "<td class=\"sidebar\" valign=\"top\">";
		admin_make_sidebar(); // includes td and /td already
		echo "</td>";
		echo "</tr>
			  </table> <!-- end super table -->
			  </div> <!-- end content -->";
		make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
				
		mysql_close($connection);
		}
	else
		{
		// error found
		
		// get title
		$title = "$config_sitenae - Error";
				
		// print html headers
		include_once('../includes/header.php');
		admin_nav();
				
		// begin HTML page content
		echo "<div id=\"content_narrow\">";
		echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
				<tr><td valign=\"top\">";

		// print as list
		echo "<img src=\"../images/warning.png\" alt=\"\">The following errors were encountered:";
		echo '<br>';
		echo '<ul>';
		for ($x=0; $x<sizeof($errorList); $x++)
			{
			echo "<li>$errorList[$x]";
			}
		echo "<form action=\"" . $_SERVER['PHP_SELF'] . "\" method=\"post\">";
		foreach ($_POST as  $key=>$value) { 
			if ($key == "intron_start" || $key == "intron_end") {
					for ($i = 1; $i <= $numIntrons; $i++){
						echo "<input type=\"hidden\" name=\"intron_start[]\" value=\"$intron_start[$i]\">" ;
						echo "<input type=\"hidden\" name=\"intron_end[]\" value=\"$intron_end[$i]\">" ;
					}
			}
			else {
				echo "<input type=\"hidden\" name=\"$key\" value=\"$value\">" ;
			}
		}
		echo "<input type=\"submit\" name=\"updateNewIntrons\" value=\"Go Back And Edit Here\" /></form>";
		echo "</ul></td>";
		echo "<td class=\"sidebar\" valign=\"top\">";
		admin_make_sidebar(); // includes td and /td already
		echo "</td>";
		echo "</tr>
			  </table> <!-- end super table -->
			  </div> <!-- end content -->";
		make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
		}
	}


// #################################################################################
// direct access - view gene table
// #################################################################################
elseif (!$_GET['new'] && !$_POST['submitNew'] && !$_POST['submitNoNew'] &&  !$_GET['geneCode'] ) {
	// get title
	$title = "$config_sitename - Gene list";
			
	// print html headers
	include_once('../includes/header.php');
	admin_nav();
			
	// begin HTML page content
	echo "<div id=\"content_narrow\">";
	echo "<table border=\"0\" width=\"850px\"> <!-- super table -->
			<tr><td valign=\"top\">";
	if( $mask_url == "true" ) {
		echo "<a href='" .$base_url . "/home.php' onclick=\"return redirect('add_gene.php?new=new');\"><b>Add gene</b></a><br />";
	}
	else {
		echo "<a href='" .$base_url . "/admin/add_gene.php?new=new'><b>";
		echo "Add gene</b></a><br />";
	}

	// print as list
	// open database connection
	@$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
	
	//select database
	mysql_select_db($db) or die ('Unable to content');
	if( function_exists(mysql_set_charset) ) {
		mysql_set_charset("utf8");
	}
	// generate and execute query from genes table
	$query = "SELECT id, geneCode, length, description, timestamp FROM ". $p_ . "genes ORDER BY geneCode";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	
	// if records present
	if (mysql_num_rows($result) > 0) {
		// iterate through result set
		// print article titles
		echo "<h1>Existing genes:</h1>\n";

		echo "<b>Create a definition for genes. Specify \"Reading frame\" if you";
		echo " want to create datasets by codon positions.</b>";

		echo "<ul>";
		while ($row = mysql_fetch_object($result)) {
			$descrutf8 = "";
			$descrutf8 = utf8_decode($row->description);

			echo "<li>";
			if( $mask_url == "true" ) {
				echo "<a href='" . $base_url . "/home.php' onclick=\"return ";
				echo "redirect('add_gene.php?geneCode=$row->geneCode');\">";
				echo "<b>$row->geneCode</b></a>";
				echo " <i>$descrutf8";
				echo ' - ' . $row->length . 'bp.';
				echo "</i>";
			}
			else {
				echo "<a href='" . $base_url . "/admin/add_gene.php?geneCode=";
				echo $row->geneCode . "'><b>$row->geneCode</b></a>";
				echo " <i>$descrutf8";
				echo ' - ' . $row->length . 'bp.';
				echo "</i>";
			}
			echo "</li>";
		}
		echo "</ul>";
	}

	// if no records present
	// display message
	else
		{
		?>
	
		<b>Create a definition for genes. Specify "Reading frame" if you want to create datasets by codon positions.</b>
		<br />
		<br />
		<font size="-1">No records currently available</font>
	
		<?php
		}
	
// close database connection
mysql_close($connection);
?>
</ul>
</td>
<?php
		echo "<td class=\"sidebar\" valign=\"top\">";
		admin_make_sidebar(); 
		echo "</td>";
		echo "</tr>
			  </table> <!-- end super table -->
			  </div> <!-- end content -->";
		make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);

	}
else
{
	{
	echo "<div id=\"rest1\"><img src=\"images/warning.png\" alt=\"\" /><span class=\"text\"> Some kind of error ocurred, but I do not know what it is, please try again!</span></div>";
	}
}
	?>

<script>
	$('#delete_gene').button({icons:{primary: 'ui-icon-alert'}}).addClass('ui-state-error');
</script>
	
</body>
</html>
