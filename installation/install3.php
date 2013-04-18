<?php
// include
require_once ('functions.php');

#print_r($_POST);exit(0);

function docheck() {
	$check_url = getParam( $_POST, 'url', '');

	$check_local_folder = getParam( $_POST, 'local_folder', '');

	$checkdatabase_host = getParam( $_POST, 'database_host', '');
	$a = clean_string($checkdatabase_host);
	$checkdatabase_host = $a[0];

	$checkdatabase_name = getParam( $_POST, 'database_name', '');
	$a = clean_string($checkdatabase_name);
	$checkdatabase_name = $a[0];

	$checkdatabase_prefix = getParam( $_POST, 'prefix', '');
	$a = clean_string($checkdatabase_prefix);
	$checkdatabase_prefix = $a[0];

	$checkdatabase_username = getParam( $_POST, 'database_username', '');
	$a = clean_string($checkdatabase_username);
	$checkdatabase_username = $a[0];

	$checkdatabase_password = $_POST['database_password'];

	$checksiteName = getParam( $_POST, 'siteName', '');
	$a = clean_string($checksiteName);
	$checksiteName = $a[0];

	if ( $checksiteName == "" ||
		 $checkdatabase_host == "" || 
		 $checkdatabase_name == "" ||
		 $checkdatabase_username == "" ||
		 $checkdatabase_password == "" ) {
		echo "<html><head><title>Error</title><link rel=\"stylesheet\" href=\"install.css\" type=\"text/css\" /></head><body><div class=\"error\"><h2><img src=\"warning.png\" alt=\"\" />You entered invalid or empty details, please try again</h2></div></body></html>";
		exit;
	}

	// try to connect to MySQL using host, user and passwd to see whether the info given is correct
	$connection = @mysql_connect($checkdatabase_host, $checkdatabase_username, $checkdatabase_password);
	if( !$connection ) {
		$error_code = "";
		$error_code = mysql_errno();
		if( $error_code != 0 ) {
			echo "<html><head><title>Error</title><link rel=\"stylesheet\" 
				  href=\"install.css\" type=\"text/css\" /></head>
				  <body>";
			echo "Error number: " . $error_code;
			echo "<div class=\"error\">
				  <h2><img src=\"error.png\" alt=\"\" /> Unable to connect to MySQL database (possible errors):</h2>";

			if( $error_code == "2002" ){
				echo "<ul>
						<li>" . mysql_error() . "</li>";
				echo "  <li>Please modify your MySQL installation file 
					      <code>/usr/local/mysql/support-files/my-large.cnf</code>
						  file parameters:
						  <ol>
							<li>Modify the lines <code>/var/mysqld/mysqld.sock</code> to <code>/tmp/mysql.sock</code></li>
							<li>Save the file as <code>/etc/my.cnf</code> and <code>/etc/mysql/my.cnf</code> in your computer</li>
							<li>Do the same to the file <code>php.ini</code>:
								<ol>
									<li><code>mysql.default_socket = /tmp/mysql.sock<code></li>
									<li>See more info here: 
									<a href='http://carlosp420.github.io/VoSeq/#installation-faq'>Installation FAQ</a>.</li>
								</ol>
							<li>Restart your computer and continue this installation (press F5 to refresh).</li>
						  </ol>
						 </li>
					</ul>";
			}
			elseif( $error_code == "1045" ){
				echo "<ul>
						<li>" . mysql_error() . "</li>";
				echo "  <li>Make sure you entered the right password for MySQL</li>
					</ul>";
			}
			else {
				echo "A unforseen error ocurred. The error number is <u>" . $error_code . "</u> 
					  <br /><br />Please report this error in: <br />
					  <a href='https://github.com/carlosp420/VoSeq/issues'>https://github.com/carlosp420/VoSeq/issues</a>";
				echo "<ul><li>" . mysql_error() . "</li></ul>";
			}


			echo "</div></body></html>";
			exit(0);
		}

	}

	mysql_query("CREATE DATABASE IF NOT EXISTS $checkdatabase_name");
	mysql_select_db($checkdatabase_name) or die ("<html><head><title>Error</title><link rel=\"stylesheet\" href=\"install.css\" type=\"text/css\" /></head><body><div class=\"error\"><h2><img src=\"error.png\" alt=\"\" /> Unable to select database</h2></div><ul><li>" . mysql_error() . "</li><li>Error number: " . mysql_errno() . "</li></ul></body></html>");


	return array("siteName" => $checksiteName,
				 "url" => $check_url,
				 "local_folder" => $check_local_folder,
				 "host" => $checkdatabase_host, 
				 "user" => $checkdatabase_username,
				 "pass" => $checkdatabase_password,
				 "prefix" => $checkdatabase_prefix,
				 "db"   => $checkdatabase_name
			);
}

$variables = docheck();
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<title>VoSeq - Web Installer</title>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<link rel="shortcut icon" href="../favicon.ico" />
<link rel="stylesheet" href="install.css" type="text/css" />

<script type="text/javascript" src="jquery.js"></script>
<script type="text/javascript" src="config.js"></script>

</head>

<body>

	
<div id="wrapper">
	<div id="header">
		<div id="joomla"><img src="header_install.png" alt="Installation" /></div>

	</div>
</div>
<div id="ctr" align="center">
	<form action="install4.php" method="post" name="form" id="form" onsubmit="return check();">
	<div class="install">
		<div id="stepbar">
			<div class="step-off">
				Installation
			</div>
			<div class="step-off">
				step 1
			</div>
			<div class="step-off">
				step 2
			</div>
			<div class="step-on">
				step 3
			</div>
			<div class="step-off">
				step 4
			</div>
			<div class="step-off">
				step 5
			</div>
		</div>
		<div id="right">
			<div class="far-right">
				<input class="button" type="submit" name="next" value="Next >>" />
  			</div>
	  		<div id="step">

	  			step 3
	  		</div>
  			<div class="clr"></div>
  			<h1>Administrator account</h1>
	  		<div class="install-form">
				<table class="content2" border="0">
				<tr>
					<td colspan="2">
						<br />
						<label class="config-label"><b>Do you want to mask URL addresses and links?:</b></label>
						<div class="mw-help-field-container">
							<span class="mw-help-field-hint">help</span>
							<span class="mw-help-field-data">
								<p>If yes, you will have URLs like this: <b>http://www.mywebsite.com/home.php</b></p>
								<p>instead of <b>http://www.mywebsite.com/story.php?code=MyVoucher001</b></p>
							</span>
						</div>
						<br />
						<input type="radio" name="mask_url" value="true"> Yes<br />
						<input type="radio" name="mask_url" value="false" checked="checked"> No
					</td>
				</tr>
				<tr>
					<td colspan="2">
						<br />
						<label class="config-label"><b>Your name:</b></label>
						<div class="mw-help-field-container">
							<span class="mw-help-field-hint">help</span>
							<span class="mw-help-field-data">
								<p>Enter your first and last names here, for example "Bob Smith".</p>
							</span>
						</div>
						<br />
						<input class="inputbox" type="text" name="admin_name" value="" size="30" />
					</td>
				</tr>
				<tr>
					<td colspan="2">
						<br />
						<label class="config-label"><b>Your login:</b></label>
						<div class="mw-help-field-container">
							<span class="mw-help-field-hint">help</span>
							<span class="mw-help-field-data">
								<p>Enter your preferred login username here, for example "bob_smith". It should not contain spaces. This is the login name you will use to log in to the system.</p>
							</span>
						</div>
						<br />
						<input class="inputbox" type="text" name="admin_login" value="" size="30" />
					</td>
				</tr>
				<tr>
					<td colspan="2">
						<br />
						<label class="config-label"><b>Password:</b></label>
						<br />
						<input class="inputbox" type="password" name="admin_password1" value="" size="30" />
					</td>
				</tr>
				<tr>
					<td colspan="2">
						<br />
						<label class="config-label"><b>Password again:</b></label>
						<br />
						<input class="inputbox" type="password" name="admin_password2" value="" size="30" />
					</td>
				</tr>
				</table>
			</div>
 
		</div>
		<div class="clr"></div>
	</div>
	<?php
		foreach($variables as $k => $v) {
			echo "<input type=\"hidden\" name=\"$k\" value=\"$v\" />\n";
		}
	?>
	</form>
</div>
<div class="clr"></div>
<div class="ctr">
	<a href="http://carlosp420.github.io/VoSeq/" target="_blank">VoSeq</a></div>

</body>
</html>
