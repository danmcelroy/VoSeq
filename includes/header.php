<?php
// #################################################################################
// #################################################################################
// Voseq includes/header.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: Builds the page header and dojo and jscript options
// #################################################################################
#setting character set in HTTP headers
header('Content-type: text/html; charset=utf8');
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
	<title>
	<?php 
		if( isset($title) ) {
			echo "$title"; 
		}
		else {
			echo "$config_sitename";
		}
	?>
	</title>
	
	<?php
	echo "<link rel=\"stylesheet\" href=\"";
	if( isset($admin) && $admin != false ) {
		echo $base_url . "/";
	}
	echo $currentTemplate . "/jquery-ui/start/jquery-ui-1.10.2.custom.css\" type=\"text/css\">";


	echo "<link rel=\"stylesheet\" href=\"";

	if( isset($admin) && $admin != false ) {
		echo $base_url . "/admin/" . $currentTemplate . "/css/1.css\" type=\"text/css\" />"; //"/admin/" 
	}
	elseif( isset($in_includes) && $in_includes != false ) {
		echo $base_url . "/" . $currentTemplate . "/css/1.css\" type=\"text/css\" />";
	}
	else {
		echo $currentTemplate . "/css/1.css\" type=\"text/css\" />";
	}

	if( isset($loginmodule) && $loginmodule != false ) {
		echo "\n<link rel=\"stylesheet\" href=\"" . $base_url . "/login/loginmodule.css\" rel=\"stylesheet\" type=\"text/css\" />";
	}

	if (isset($yahoo_map) && $yahoo_map != false) {
		echo "\n<script src=\"http://api.maps.yahoo.com/ajaxymap?v=3.7&appid=" . $yahoo_key . "\">";
		echo "</script>";
	}
	echo "\n\n<script type=\"text/javascript\" src=\"" . $base_url . "/includes/jquery.js\"></script>\n";
	echo "\n\n<script type=\"text/javascript\" src=\"" . $base_url . "/includes/jquery-ui.js\"></script>\n";
	echo "\n\n<script type=\"text/javascript\" src=\"" . $base_url . "/includes/jquery_download.js\"></script>\n";

    if( isset($plupload) ) {
		echo "\n<link rel=\"stylesheet\" href=\"plupload/js/jquery.plupload.queue/css/jquery.plupload.queue.css\" rel=\"stylesheet\" type=\"text/css\" />";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/plupload.full.js\"></script>\n";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/jquery.plupload.queue/jquery.plupload.queue.js\"></script>\n";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/plupload.js\"></script>\n";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/plupload.gears.js\"></script>\n";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/plupload.silverlight.js\"></script>\n";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/plupload.flash.js\"></script>\n";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/plupload.browserplus.js\"></script>\n";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/plupload.html4.js\"></script>\n";
	    echo "<script type=\"text/javascript\" src=\"" . $base_url . "/admin/plupload/js/plupload.html5.js\"></script>\n";
    }
	?>
	
	<link rel="SHORTCUT ICON" href="<?php echo $base_url . "/favicon.ico"; ?>" />
	<meta content="Gvim" name="GENERATOR" />
	<meta content="Carlos Pe&ntilde;a" name="author" />
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />

<?php 
if( isset($dojo) && $dojo == true ) {
	echo "<script type=\"text/javascript\" src=\"";
	if ($admin) {
		echo "../";
		}
	echo "dojo/dojo.js\"></script>\n<script type=\"text/javascript\">\n";
	echo "\tdojo.require(\"dojo.widget.*\");\n"; // load code relating to widget managing fuctions
	foreach ($whichDojo as $value) {
		echo "\tdojo.require(\"dojo.widget." . $value . "\");\n"; // put_dojo();
		}
	echo "</script>\n";
}

// redirect function for masking urls
if ($mask_url == "true") {
	echo "
	<script type=\"text/javascript\">
		function redirect(URL) {
			document.location=URL;
			return false;
		}
	</script>\n";
}
	
if ( !isset($delete_buttons) ) {
echo "
<script type=\"text/javascript\">
	$(document).ready(function() {
		$('#delete_record').click(function(){
			var id = $('#delete_record').attr('name');
			if(confirm('Delete record? Notice that you will lose also accompanying sequences!')) {
				$.post('delete_record.php', {id: id},
					function(data) {
						if(data == 'ok') {
							alert('Your record was successfully deleted');
							window.location.replace('admin.php');
						}
						else {
							alert('I could not remove your record');
						}
					});
			}
		});


		$('#delete_gene').click(function(){
			var id = $('#delete_gene').attr('name');
			if(confirm('Delete gene? Notice that you will lose also accompanying sequences!')) {
				$.post('delete_gene.php', {id: id},
					function(data) {
						if(data == 'ok') {
							alert('Your gene was successfully deleted');
							window.location.replace('admin.php');
						}
						else {
							alert('I could not remove your gene');
						}
					});
			}
		});
   

		$('#delete_sequence').click(function(){
			var id = $('#delete_sequence').attr('name');
			if(confirm('Delete sequence?')) {
				$.post('delete_sequence.php', {id: id},
					function(data) {
						if(data == 'ok') {
							alert('Your sequence was successfully deleted');
							window.location.replace('admin.php');
						}
						else {
							alert('I could not remove your sequence');
						}
					});
			}
		});
   


    $('a.delete').on('click',function(e){
        e.preventDefault();
		if( confirm('Delete photo?') ) {
			var imageID = $(this).parent('.voucher')[0].id;
			var voucherImage = $(this).siblings().attr('href');
			var thumbnail = $(this).siblings().children().attr('src');
			$.post('delete_photo.php', { 'voucherImage': voucherImage,
									 'thumbnail': thumbnail
									 });
			$(this).closest('.voucher')
				.fadeTo(300,0,function(){
					$(this)
						.animate({width:0},200,function(){
							$(this).remove();
						});
				});
			}
		});
	});
</script>\n";
}


?>
</head>
<body>
