<?php
// ############################################################################
// ############################################################################
// Voseq admin/admin.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: splits a big file in chuncks to be uploaded by upload.php
//                  using plupload
//
// ############################################################################



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

$plupload = true;

$delete_buttons = false;

$title = $config_sitename . "| Database merger";
  
include_once('../includes/header.php');

admin_nav();

// scan upload dir and remove files
$scanfiles = scandir('uploads');
foreach ($scanfiles as $num => $scanfile){
	if ($scanfile !="." && $scanfile !=".."){
		unlink("uploads/".$scanfile);
	}
}


echo "<div id=\"content_narrow\" style='font-size: 12px;'>";
echo "<h1>Merge your database with another SQL file</h1>";
echo "<p>Upload another MySQL database file that were generated either by:
    <ul>
        <li>Using VoSeq's <b>\"Backup database\"</b> button.</li>
        <li>From a terminal or console with the command: </li>
        <ul>
            <li><code>mysqldump " . $db . " -uroot -pmy_password > 
                db-file.sql</code></li>
        </ul>
        <li>Note that <b><u>no data will be replaced, but new data will be added</u></b> from the 
             uploaded file.</li>
    </ul>
    </p>";


$form_output = "
<div id='container'>
    <div id='filelist'>No runtime found. <b>Your web browser doesn't support uploading chunked files!</b>
    Please use a modern browser such as Firefox or Google Chrome to upload files.
    </div>
    <a class='plupload_button plupload_add' id='pickfiles' href='javascript:;'>Select file</a>
    <a class='plupload_button plupload_start' id='uploadfiles' href='javascript:;'>Upload file</a>
</div>
";


echo $form_output;
echo "</div>";
make_footer($date_timezone, $config_sitename, $version, $base_url, $p_);
?>

<script type="text/javascript">
function mysqlmerge(filename) {
	var r=confirm("File uploaded.\n\nThis merger will add new data\nfrom the chosen db.sql to your\nexisting database!\n\nBe sure to check the result output after!\n\nContinue processing?");
	if (r==false){x="Aborted";}
	else {
	  //x="You pressed Cancel!";
		jQuery.noConflict();
		jQuery.post('mysqlmerger.php', { name: filename })
		//$.download('mysqlmerger.php', 'filename='+filename, '');
			//.done(function(data) {html(data);});
			 .done(function(data) { 
					alert("Finished processing!\nCheck and copy output from main screen!!");
					// var r=confirm("-----------------------\n" + data.substring(0,1000) + "\n\n...\n-----------------------\n" + "Save merge result output to file?");
					  // if (r==true){
							 $('filelist').innerHTML = "<div><br>"+ data +"<br><br></div>";
						 // }
					  // else{}
			 })
			.fail(function() { alert("Error in processing. Try again or contact VoSeq creators."); })
	 }
}

// Custom example logic
function $(id) {
	return document.getElementById(id);	
}

var uploader = new plupload.Uploader({
	runtimes : 'gears,html5,flash,silverlight,browserplus',
	browse_button : 'pickfiles',
	container: 'container',
    chunk_size: '1mb',
	max_file_size : '100mb',
	url : 'upload.php',
	resize : {width : 320, height : 240, quality : 90},
	flash_swf_url : 'plupload/js/plupload.flash.swf',
	silverlight_xap_url : 'plupload/js/plupload.silverlight.xap',
	filters : [
		{title : "SQL files", extensions : "sql"}
	]
});

uploader.bind('Init', function(up, params) {
	$('filelist').innerHTML = "<div>Your web browser allows uploading big files using " + params.runtime + "!</div>";
});

uploader.bind('FilesAdded', function(up, files) {
	for (var i in files) {
		$('filelist').innerHTML += '<div  id="' + files[i].id + '">' + files[i].name + ' (' + plupload.formatSize(files[i].size) + ') <b></b></div>';
	}
});

uploader.bind('UploadProgress', function(up, file) {
	$(file.id).getElementsByTagName('b')[0].innerHTML = '<span id="progress">' + file.percent + "%</span>";
});

// my code
uploader.bind('FileUploaded', function(up, file, resp) {
    var my_json = jQuery.parseJSON(resp.response);
    if( my_json.result == null ) {
		$('filelist').innerHTML = "<div STYLE=\"font-size: 16px;\"><br><br><b>Wait for it... ... ...<b><br><br></div>";
        mysqlmerge(file.name);
    }
});

$('uploadfiles').onclick = function() {
	uploader.start();
	return false;
};

uploader.init();

</script>


</body>
</html>
