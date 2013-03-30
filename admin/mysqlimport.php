<?php
// ############################################################################
// ############################################################################
// Voseq admin/admin.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Entry page for administrator interface
//
// ############################################################################


/*
 * taken from http://dan.cx/blog/2006/12/restore-mysql-dump-using-php
 *
 * Restore MySQL dump using PHP
 * (c) 2006 Daniel15
 * Last Update: 9th December 2006
 * Version: 0.2
 *
 * modified Carlos Peña 2013-03-29
 */


//check admin login session
include'../login/auth-admin.php';
// includes
ob_start();//Hook output buffer - disallows web printing of file info...
include'../conf.php';
ob_end_clean();//Clear output buffer//includes
include'../functions.php';
include'adfunctions.php';
include'admarkup-functions.php';
#include'../login/redirect.html';

// to indicate this is an administrator page
$admin = true;

$dojo = false;

$title = $config_sitename . " database importer";
  
include_once('../includes/header.php');

admin_nav();

echo "<script>";
echo "$(document).ready(function() {
    $('#droppable').droppable({
            activeClass: 'ui-state-hover',
            hoverClass: 'ui-state-active',
            drop: function( event, ui) {
                $(this)
                    .addClass('ui-state-highlight')
                    .find('span')
                    .html('Dropped!');
                }
        });
    });";
echo "</script>";

echo "<div id=\"content_narrow\">";
echo "<h1>Restore your database from a backup SQL file</h1>";

echo "<div id='droppable' class='ui-widget-header'>";
echo "<p>accept: </p>";
echo "<button id='upload_droppable'>";
echo "<img class='upload_icon' src='images/pixel.gif'>";
echo "<span id='upload-button-text'>Select file to upload</span>";
echo "<span id='drag-and-drop'>";
echo "or drag a file from your desktop and drop it here</span>";
echo "</button>";
echo "</div>";

echo "</div>";
make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
echo "</body>\n</html>";
exit(0);





// Name of the file
$filename = 'test.sql';
 
 
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
 
?>
