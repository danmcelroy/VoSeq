<?php

// Name of the file
$filename = "uploads/" . trim($_POST['name']);
$filename = preg_replace('/[^\w\._^\/]+/', '_', $filename);
if( !file_exists($filename) ) {
    echo "Error, couldn't find file: ".$filename."!";
    exit(0);
}

//check admin login session
include'../login/auth-admin.php';
// includes
ob_start();//Hook output buffer - disallows web printing of file info...
include'../conf.php';
ob_end_clean();//Clear output buffer//includes

// Connect to MySQL server
mysql_connect($host, $user, $pass) or die('Error connecting to MySQL server: ' . mysql_error());
// Select database
mysql_select_db($db) or die('Error selecting MySQL database: ' . mysql_error());
 
// Temporary variable, used to store current query
$templine = '';
// array for collection of table names
$tablenames = array();
// Read in entire file
$lines = file($filename);
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
		if (strpos($tn, "members") !== FALSE) {$line = str_replace($tn, $p_ . "members", $line); }
		if (strpos($tn, "vouchers") !== FALSE) {$line = str_replace($tn, $p_ . "vouchers", $line); }
		if (strpos($tn, "genesets") !== FALSE) {$line = str_replace($tn, $p_ . "genesets", $line); }
		elseif (strpos($tn, "genes") !== FALSE && strpos($tn, "genesets") === FALSE) {$line = str_replace($tn, $p_ . "genes", $line); }
		if (strpos($tn, "sequences") !== FALSE) {$line = str_replace($tn, $p_ . "sequences", $line); }
		if (strpos($tn, "primers") !== FALSE) {$line = str_replace($tn, $p_ . "primers", $line); }
		if (strpos($tn, "taxonsets") !== FALSE) {$line = str_replace($tn, $p_ . "taxonsets", $line); }
		if (strpos($tn, "search_results") !== FALSE) {$line = str_replace($tn, $p_ . "search_results", $line); }
		elseif (strpos($tn, "search") !== FALSE && strpos($tn, "search_results") === FALSE) {$line = str_replace($tn, $p_ . "search", $line); }
		//if ( strpos($line,'CREATE TABLE') !== FALSE ){echo substr($line,0,50) . "\n";}
		}
    // Add this line to the current segment
    $templine .= $line;
    // If it has a semicolon at the end, it's the end of the query
    if (substr(trim($line), -1, 1) == ';')
    {
        // Perform the query
        mysql_query($templine) or print('Error performing query \'<strong>' . $templine . '\': ' . mysql_error() . '\n');
        // Reset temp variable to empty
        $templine = '';
    }
}
 
if( file_exists($filename) ) {
    unlink($filename);
    echo "Finished importing database!";
}

?>
