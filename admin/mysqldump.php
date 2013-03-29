<?php

session_start();
//Check whether the session variable SESS_MEMBER_ID is present or not
if(!isset($_SESSION['SESS_MEMBER_ID']) || (trim($_SESSION['SESS_MEMBER_ID']) == '')) {
	header("location: ../home.php");
	exit(0);
}
elseif( isset($_SESSION['SESS_ADMIN']) && $_SESSION['SESS_ADMIN'] == '0') {
	header("location: ../login/admin-failed.php");
	exit(0);
}


//includes
ob_start();//Hook output buffer - disallows web printing of file info...
include '../conf.php';
if( $mask_url == "true" ) {
	ob_end_clean();//Clear output buffer
}
else {
	ob_clean();
}

error_reporting(E_ALL);

/* 
 * try to get the base directory for mysql
 */
$mysql_dir = "";
@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect'. mysql_error());
mysql_select_db('information_schema');
$query = "select variable_value from global_variables where variable_name='BASEDIR'";
$result = mysql_query($query);
if( mysql_num_rows($result) > 0 ) {
	while( $row = mysql_fetch_row($result) ) {
		$mysql_dir = $row[0];
		break;
	}
}


if( $mysql_dir == "" ) {
	ob_start();
	phpinfo();
	$s = ob_get_contents();
	ob_end_clean();
	
	$s = explode("\n", $s);
	
	foreach( $s as $line ) {
		preg_match("/MYSQL_INCLUDE =>(.+)/", $line, $match);
		if( $match ) {
			$mysql_dir = trim($match[1]);
			$mysql_dir = preg_replace("/^-I/", "", $mysql_dir);
			$mysql_dir = preg_replace("/\/+include$/", "", $mysql_dir);
		}
	}
}




/* 
 * now do the mysqldump!
 */
if( $mysql_dir == "" ) {
	// hope that mysqldump is in the path
	$cmd = "mysqldump";
}
else {
	$cmd = $mysql_dir . "/bin/mysqldump";
}

if( strtoupper(substr(PHP_OS, 0, 3)) === "WIN" ) {
	$cmd .= '.exe"';
	$cmd = str_replace("/", "\\", $cmd);
	$cmd = preg_replace('/^/', '"', $cmd);
}

@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect'. mysql_error());
$cmd .= " " . $db;
$cmd .= " -h" . $host;
$cmd .= " -u" . $user;
$cmd .= " -p" . $pass;
	
//$f = fopen("a.txt", "w");
//fwrite($f, $cmd); exit(0);

unset($output);
exec($cmd, $output);

date_default_timezone_set($date_timezone);

//save file
$filename  = 'db-backup-';
$filename .= date('Y-m-d') . "_" ;
$filename .= time() . ".sql";

$dump = "";
foreach( $output as $item) {
	$dump .= $item . "\n";
}

header('Content-Description: File Transfer');
header('Content-Type: text/plain');
header('Content-Disposition: attachment; filename=' . basename($filename));
header('Expires: 0');
header('Pragma: public');
echo $dump;
?>
