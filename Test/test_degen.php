<?php
require_once(dirname(__FILE__) . '/simpletest/autorun.php');
require_once(dirname(__FILE__) . '/../includes/translation_functions.php');

class TestDegen extends UnitTestCase {
    function test_degen_coding() {
        # $seq, $genetic_code, $rframe
        $seq         = "aATGATCTACAAaTGCGGTGGTATCGATAAACGTACCATTGAGAAATTCGAGAAGGAG";
        $encoded_seq = "AATGATHTAYAARTGYGGNGGNATHGAYAARMGNACNATHGARAARTTYGARAARGAR";
        $genetic_code = "1";
        $rframe = "2";
        $this->assertEqual(degen_coding($seq, $genetic_code, $rframe), $encoded_seq);
    }
}

?>
