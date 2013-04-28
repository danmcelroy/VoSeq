<?php
// #################################################################################
// #################################################################################
// Voseq login/account_change.php
// Copyright (c) 2006, PHPSense.com
// All rights reserved.
//
// Modified by Carlos Peña & Tobias Malm
// Mods: added admin-user choice and include() output buffer
//
// Script overview: prints the account change form
// #################################################################################

//header('Content-type: text/html; charset=utf8');

//check admin login session
include'auth.php';
//includes

ob_start();//Hook output buffer - disallows web printing of file info...
include'../conf.php';
ob_end_clean();//Clear output buffer//includes
include'../functions.php';
//include'../admin/adfunctions.php';
include'../markup-functions.php';
//include'../admin/admarkup-functions.php';
error_reporting (E_ALL ^ E_NOTICE);

// admin?
//$admin = true;

// process title
$title = $config_sitename;

// loginmodule
$loginmodule = true;

// print html headers
include_once'../includes/header.php';

// print navegation bar
nav();

// header: send bugs to me message
standardHeader($title, $intro_msg);


if( isset($_SESSION['ERRMSG_ARR']) && is_array($_SESSION['ERRMSG_ARR']) && count($_SESSION['ERRMSG_ARR']) >0 ) {
	echo '<ul  class="err">';
	foreach($_SESSION['ERRMSG_ARR'] as $msg) {
		echo '<li style="text-align:center">',$msg,'</li>'; 
	}
	echo '</ul>';
	unset($_SESSION['ERRMSG_ARR']);
}

// open database connections
@$connection = mysql_connect($host, $user, $pass) or die('Unable to connect');
mysql_select_db($db) or die ('Unable to select database');
if( function_exists(mysql_set_charset) ) {
	mysql_set_charset("utf8");
}
// generate and execute query from members table
$sessid = $_SESSION['SESS_MEMBER_ID'];
$queryac = "SELECT member_id, firstname, lastname, login FROM " . $p_ . "members WHERE member_id=$sessid";
$resultac = mysql_query($queryac) or die("Error in query: $query. " . mysql_error());
// if records present
if (mysql_num_rows($resultac) > 0) {
	while ($row = mysql_fetch_object($resultac)) {
	if (isset($row->firstname)) {$firstn = $row->firstname;} else {$firstn = '';}
	if (isset($row->lastname)) {$lastn = $row->lastname;} else {$lastn = '';}
	if (isset($row->login)) {$logn = $row->login;} else {$logn = '';}
	}
}

?>
<form id="loginForm" name="loginForm" method="post" action="account_change_exec.php" accept-charset="UTF-8">
  <table width="300" border="0" align="center" cellpadding="2" cellspacing="0">
	<tr><td colspan=2><h1><b>Change account information/password</b></h1></td></tr>
	<tr>
	  <th>First Name </th>
	  <td><input name="fname" type="text" class="textfield" id="fname" value="<?php echo $firstn; ?>"/></td>
	</tr>
	<tr>
	  <th>Last Name </th>
	  <td><input name="lname" type="text" class="textfield" id="lname" value="<?php echo $lastn; ?>"/></td>
	</tr>
	<tr>
	  <th width="124">Login</th>
	  <td width="168"><input name="login" type="text" class="textfield" id="login" value="<?php echo $logn; ?>"/></td>
	</tr>
	<tr>
	  <th>Old password</th>
	  <td><input name="opassword" type="password" class="textfield" id="opassword" /></td>
	</tr>
	<tr>
	  <th>Password</th>
	  <td><input name="password" type="password" class="textfield" id="password" /></td>
	</tr>
	<tr>
	  <th>Confirm Password </th>
	  <td><input name="cpassword" type="password" class="textfield" id="cpassword" /></td>
	</tr>
	<br />
	<tr><input type="hidden" name="memberid" value="<?php echo $taxonadds2; ?>" >
		<td>&nbsp;</td>
		<td><input type="submit" name="Submit" value="Update user account" /></td>
	</tr>
</table>
</form>
</div> <!-- end content -->

<!-- standard page footer begins -->
<?php
make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
?>

</body>
</html>
