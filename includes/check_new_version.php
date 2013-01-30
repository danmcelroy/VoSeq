<?php

include("../functions.php");

// #################################################################################
// Section: index.php
//
// Using github API v3 to get all tags in repository
// @output: most recent tag. Example: 1.3.7
//          if failure, returns empty ""
// #################################################################################
function check_repo_tags() {
	$url = "https://api.github.com/repos/carlosp420/VoSeq/git/refs/tags";

	$tags = get_from_URL($url);
	$tags = json_decode($tags);

	if( $tags != NULL ) {
		$most_recent = array_pop($tags);
		$most_recent = $most_recent->ref;
		preg_match("/refs\/tags\/v(.+)/", $most_recent, $match);
	
		$most_recent = $match[1];
		return $most_recent;
	}
	else {
		return "";
	}
}

?>
