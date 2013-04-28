<?php
// #################################################################################
// #################################################################################
// Voseq login/register-success.php
// Copyright (c) 2006, PHPSense.com
// All rights reserved.
//
// Modified by Carlos Peña & Tobias Malm
// Mods: inserted include() output buffer and admin-link
//
// Script overview: Print account change success page
// #################################################################################
//check login session
include'auth.php';
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-9" />
<title>Account change successful</title>
<link href="loginmodule.css" rel="stylesheet" type="text/css" />
</head>
<body>
<h1>Account change successful</h1>
<?php
include '../conf.php';
if( $mask_url == "true" ) {
	include 'redirect.html';
	echo "<p><a href='" . $base_url . "/home.html'  onclick=\"return redirect('../index.php');\" >Click here</a> to continue your work.</p>";
}
else {
	echo "<p><a href='" . $base_url . "/index.php';\" >Click here</a> to continue your work.</p>";
}
?>
</body>
</html>
