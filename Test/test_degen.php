<?php
require_once(dirname(__FILE__) . '/simpletest/autorun.php');
require_once(dirname(__FILE__) . '/../includes/translation_functions.php');

class TestDegen extends UnitTestCase {
    function test_degen_coding_reading_frame_2() {
        # $seq, $genetic_code, $rframe
        $seq         = "aATGATCTACAAaTGCGGTGGTATCGATAAACGTACCATTGAGAAATTCGAGAAGGAG";
        $encoded_seq = "AATGATHTAYAARTGYGGNGGNATHGAYAARMGNACNATHGARAARTTYGARAARGAR";
        $genetic_code = "1";
        $rframe = "2";
        $this->assertEqual(degen_coding($seq, $genetic_code, $rframe), $encoded_seq);
    }

    function test_degen_coding_reading_frame_1() {
        # $seq, $genetic_code, $rfs
        $seq         = "ATGATCTACAAaTGCGGTGGTATCGATAAACGTACCATTGAGAAATTCGAGAAGGAG";
        $encoded_seq = "ATGATHTAYAARTGYGGNGGNATHGAYAARMGNACNATHGARAARTTYGARAARGAR";
        $genetic_code = "1";
        $rframe = "1";
        $this->assertEqual(degen_coding($seq, $genetic_code, $rframe), $encoded_seq);
    }
}

?>
