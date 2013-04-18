<?php
// #################################################################################
// #################################################################################
// Voseq markup-functions.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Functions for outoput formats
// #################################################################################
// #################################################################################
// Section: Startup/includes
// #################################################################################
#include 'login/redirect.html';
ob_start();//Hook output buffer - disallows web printing of file info...
include 'conf.php';
ob_end_clean();//Clear output buffer//includes

// #################################################################################
// Section: nav() function
// This function draws the top navigation bar HOME, SEARCH, ADMIN
// #################################################################################
function nav() {
	ob_start();//Hook output buffer - disallows web printing of file info...
	include 'conf.php';
	ob_end_clean();//Clear output buffer//includes
	echo "<div id=\"menu\">";

	// masking URLs, this variable is set to "true" or "false" in conf.php file
	if( $mask_url == "true" ) {
		echo "<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/index.php');\" title='This link takes you back to the homepage' >home</a>";
		echo "<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/search.php');\" title='Go to search page' >search</a>";
		echo "<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/admin/admin.php');\" title='Go to Administration page' >admin</a>";
		echo "</div>\n\n";
	}
	else {
		echo "<a href='" . $base_url . "/index.php' title='This link takes you back to the homepage' >home</a>";
		echo "<a href='" . $base_url . "/search.php' title='Go to search page' >search</a>";
		echo "<a href='" . $base_url . "/admin/admin.php' title='Go to Administration page' >admin</a>";
		echo "</div>\n\n";
	}
}
		
// #################################################################################
// Section: standardheader() function
// This function draws standardheader and intro text
// #################################################################################
function standardHeader($title) {
echo "<!--standard page header begins-->
		<div id=\"header\">
			<h1>". $title . "</h1>
			<h2>This is <b>VoSeq</b>. VoSeq is a database to store voucher and sequence data. 
			    Please send all bug complaints to <br />Carlos Peña (<i>mycalesis@gmail.com</i>) or<br />
				Tobias Malm (<i>tobias.malm@uef.fi</i>)&nbsp;&nbsp;</h2>
			
			<p class=\"introduction\"><b>1.</b> This is the <b>user</b> interface of this database, which means you can only search and look for records, retrieve sequences, get details of vouchers and look at vouchers' pictures.
			<br /><b>2.</b> If you want to add, delete or update records and sequences, please click the <b>\"ADMIN\"</b> button (for administration).
			<br />
				<b>3.</b> The full documentation is here: <a href='http://carlosp420.github.io/VoSeq/'>http://carlosp420.github.io/VoSeq/</a>.
			<br />
				<b>4.</b> The source code is available at: <a href='https://github.com/carlosp420/VoSeq'>https://github.com/carlosp420/VoSeq</a>.
			<br />
				<b>5.</b> If you find VoSeq useful you might want cite us:
				<div class='highlight'><b>Peña, C. &amp; Malm, T.</b> 2012. VoSeq: a Voucher and DNA Sequence Web Application. <i>PLoS ONE</i>, 7(6): e39071. <a href=\"http://dx.doi.org/10.1371/journal.pone.0039071\">doi:10.1371/journal.pone.0039071</a></div>
			
			</p>
		</div>
		<!-- standard page header ends -->";
}

function end_divs() {
echo "</div> <!-- end mainbar -->
		</div> <!-- end content -->";
}

// #################################################################################
// Section: make_sidebar() function
// This functions draws the sidebar on the right with logo and link to tools
// #################################################################################
function make_sidebar() {
	ob_start();//Hook output buffer - disallows web printing of file info...
	include 'conf.php';
	ob_end_clean();//Clear output buffer//includes
	echo "<img width=\"160px\" height=\"79px\" src=\"" . $base_url . "/images/logo-small.jpg\" alt=\"VoSeq database\" class=\"logo\" />";
	echo "<h1>Tools:</h1>
			<div class=\"submenu\">";
	if( $mask_url == "true" ) {
		echo   "<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/blast_new.php');\">Blast new sequence</a>
				<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/view_table.php');\">Interactive overview table</a>
				<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/genes.php');\">View genes</a>
				<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/create_dataset.php');\">Create new dataset</a>
				<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/create_table.php');\">Create Voucher table</a>
				<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/create_gene_table.php');\">Create Gene table</a>
				<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/create_genbank_fasta_file.php');\">Create GenBank FASTA file</a>
				<a href='" . $base_url . "/home.php'  onclick=\"return redirect('". $base_url . "/share_data_gbif.php');\">Share data with GBIF</a>
				";
	}
	else {
		echo   "<a href='" . $base_url . "/blast_new.php'\">Blast new sequence</a>
				<a href='" . $base_url . "/view_table.php'\">Interactive overview table</a>
				<a href='" . $base_url . "/genes.php'\">View genes</a>
				<a href='" . $base_url . "/create_dataset.php'\">Create new dataset</a>
				<a href='" . $base_url . "/create_table.php'\">Create Voucher table</a>
				<a href='" . $base_url . "/create_gene_table.php'\">Create Gene table</a>
				<a href='" . $base_url . "/create_genbank_fasta_file.php'\">Create GenBank FASTA file</a>
				<a href='" . $base_url . "/share_data_gbif.php'\">Share data with GBIF</a>
				";
	}
			
	echo "</div>";
}



// #################################################################################
// Section: make_footer() function
// This functions draws the footer
// #################################################################################
function make_footer($date_timezone, $config_sitename, $version, $base_url) {
	ob_start();//Hook output buffer - disallows web printing of file info...
	include 'conf.php';
	ob_end_clean();//Clear output buffer//includes
			
	//fixing some output variables
	// open database connection
	@$connection = mysql_connect($host, $user, $pass) or die ('Unable to connect!');
	//select database
	mysql_select_db($db) or die ('Unable to content');
	if( function_exists(mysql_set_charset)) {
		mysql_set_charset("utf8");
	}

	$num_rows_vouchers = array();
	$query = "SELECT * FROM " . $p_ . "vouchers";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	$num_rows_vouchers['all'] = mysql_num_rows($result);
	$count_array = array("orden", "family","genus", "genus, species");
	foreach ($count_array as $count){
		$query = "SELECT $count FROM ". $p_ . "vouchers GROUP BY $count";
		$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
		$num_rows_vouchers[$count] = mysql_num_rows($result);
	}
	$query = "SELECT * FROM " . $p_ . "sequences";
	$result = mysql_query($query) or die("Error in query: $query. " . mysql_error());
	$num_rows_sequences = mysql_num_rows($result);

	if( $php_version == "5" ) {
		// date_default_timezone_set($date_timezone); php5
		date_default_timezone_set($date_timezone);
	}
	
	echo "<!-- standard page footer begins -->\n<div id=\"footer\">" . date('Y') . ' ' . $config_sitename;
	echo "\n <a href=\"" . $base_url . "/home.php\" onclick=\"return redirect('" . $base_url . "/changelog.txt');\" title='check verion history' >version " . $version . "</a>";
	echo " \n Logged in as: " . $_SESSION['SESS_FIRST_NAME'] ." ". $_SESSION['SESS_LAST_NAME'] . "\n <a href='
	" . $base_url . "/home.php' onclick=\"return redirect('". $base_url . "/login/logout.php');\">logout </a>";
	echo "<br />Now with <b>". $num_rows_vouchers['all'] ."</b> vouchers, over ".$num_rows_vouchers['orden']." orders, 
	".$num_rows_vouchers['family']." familes, ".$num_rows_vouchers['genus']." genera and ".$num_rows_vouchers["genus, species"]." species,
	with together ". $num_rows_sequences ." sequences! <br />\n<img width=\"80px\" height=\"15px\" src=\"" . $base_url . "/images/colofon_xhtml.png\" alt=\"Valid XHTML\" title=\"Valid XHTML\" />
	<img width=\"80px\" height=\"15px\" src=\"" . $base_url . "/images/colofon_css.png\" alt=\"Valid CSS\" title=\"Valid CSS\" />\n</div>";
}

?>
