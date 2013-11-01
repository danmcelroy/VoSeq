<?php

if( file_exists("changelog.md") ) {
	$output = file_get_contents("changelog.md");
}

header("Content-Type: text/plain; charset=utf-8");
echo $output;
?>
