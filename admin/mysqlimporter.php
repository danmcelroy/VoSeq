<?php

// Name of the file
$filename = "uploads/" . trim($_POST['name']);
$filename = str_replace("-", "_", $filename);

if( !file_exists($filename) ) {
    echo "Error, couldn't find that file!";
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
// Read in entire file
$lines = file($filename);
// Loop through each line
foreach ($lines as $line)
{
    // Skip it if it's a comment
    if (substr($line, 0, 2) == '--' || $line == '' || substr($line, 0, 3) == '/*!') {
        continue;
    }
 
    // Add this line to the current segment
    $templine .= $line;
    // If it has a semicolon at the end, it's the end of the query
    if (substr(trim($line), -1, 1) == ';')
    {
        // Perform the query
        mysql_query($templine) or print('Error performing query \'<strong>' . $templine . '\': ' . mysql_error() . '<br /><br />');
        // Reset temp variable to empty
        $templine = '';
    }
}
 
if( file_exists($filename) ) {
    unlink($filename);
    echo "Finished importing database!";
}

?>
