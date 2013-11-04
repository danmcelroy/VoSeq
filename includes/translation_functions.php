<?php
// #################################################################################
// #################################################################################
// Voseq includes/translation_functions.php
// author(s): Carlos Peña & Tobias Malm
// license   GNU GPL v2
// source code available at https://github.com/carlosp420/VoSeq
//
// Script overview: This is a modified script retrieveing the functions from a script
//  downloaded from http://www.biophp.org/minitools/dna_to_protein/
//
// original script information:
// author    Joseba Bikandi
// license   GNU GPL v2
// source code available at  biophp.org
// #################################################################################
// #################################################################################
// Section: find_ORF() function
// #################################################################################
function find_ORF($frame, $protsize,$only_coding,$trimmed){
        foreach ($frame as $n => $peptide_sequence){
                $peptide_sequence=strtolower($peptide_sequence);
                $oligo=preg_split('/\*/',$peptide_sequence);
                foreach ($oligo as $m => $val){
                        if (strlen($val)>=$protsize){
                                if ($trimmed==1){
                                        $oligo[$m]=substr($val,0,strpos($val,"m")).strtoupper(substr($val,strpos($val,"m")));
                                }else{
                                        $oligo[$m]=strtoupper($val);
                                }
                        }
                }
                $new_peptide_sequence="";
                foreach ($oligo as $m => $val){if($m!=0){$new_peptide_sequence.="*".$val;}else{$new_peptide_sequence.=$val;}}
                // to avoid showing no coding, remove them from output sequence
                if($only_coding==1){$new_peptide_sequence=preg_replace("/f|l|i|m|v|s|p|t|a|y|h|q|n|k|d|e|c|w|r|g|x]/","_",$new_peptide_sequence);}
                $frame[$n]=$new_peptide_sequence;
        }
        return $frame;
}

// #################################################################################
// Section: RevComp_DNA() function
// #################################################################################
function RevComp_DNA($seq){
        $seq= strtoupper($seq);
        $original=  array("(A)","(T)","(G)","(C)","(Y)","(R)","(W)","(S)","(K)","(M)","(D)","(V)","(H)","(B)");
        $complement=array("t","a","c","g","r","y","w","s","m","k","h","b","d","v");
        $seq = preg_replace ($original, $complement, $seq);
        $seq= strtoupper ($seq);
        return $seq;
}
// #################################################################################
// Section: translate_DNA_to_protein() function
// #################################################################################
function translate_DNA_to_protein($seq,$genetic_code){

        // $aminoacids is the array of aminoacids
        $aminoacids=array("F","L","I","M","V","S","P","T","A","Y","*","H","Q","N","K","D","E","C","W","R","G","X");

        // $triplets is the array containning the genetic codes
        // Info has been extracted from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode

        // Standard genetic code
        $triplets[1]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG |TGA )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Vertebrate Mitochondrial
        $triplets[2]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG |AGA |AGG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG |TGA )","(CG. )","(GG. )","(\S\S\S )");
        // Yeast Mitochondrial
        $triplets[3]=array("(TTT |TTC )","(TTA |TTG )","(ATT |ATC )","(ATG |ATA )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. |CT. )","(GC. )","(TAT |TAC )","(TAA |TAG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG |TGA )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Mold, Protozoan and Coelenterate Mitochondrial. Mycoplasma, Spiroplasma
        $triplets[4]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG |TGA )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Invertebrate Mitochondrial
        $triplets[5]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC )","(ATG |ATA )","(GT. )","(TC. |AG. )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG |TGA )","(CG. )","(GG. )","(\S\S\S )");
        // Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
        $triplets[6]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TGA )","(CAT |CAC )",
                        "(CAA |CAG |TAA |TAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Echinoderm Mitochondrial
        $triplets[9]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AG. )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAA |AAT |AAC )","(AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG |TGA )","(CG. )","(GG. )","(\S\S\S )");
        // Euplotid Nuclear
        $triplets[10]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC |TGA )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Bacterial and Plant Plastid
        $triplets[11]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG |TGA )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Alternative Yeast Nuclear
        $triplets[12]=array("(TTT |TTC )","(TTA |TTG |CTA |CTT |CTC )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC |CTG )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG |TGA )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Ascidian Mitochondrial
        $triplets[13]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC )","(ATG |ATA )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG |TGA )","(CG. )","(GG. |AGA |AGG )","(\S\S\S )");
        // Flatworm Mitochondrial
        $triplets[14]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AG. )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC |TAA )","(TAG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC |AAA )","(AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG |TGA )","(CG. )","(GG. )","(\S\S\S )");
        // Blepharisma Macronuclear
        $triplets[15]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TGA )","(CAT |CAC )",
                        "(CAA |CAG |TAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Chlorophycean Mitochondrial
        $triplets[16]=array("(TTT |TTC )","(TTA |TTG |CT. |TAG )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TGA )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Trematode Mitochondrial
        $triplets[21]=array("(TTT |TTC )","(TTA |TTG |CT. )","(ATT |ATC )","(ATG |ATA )","(GT. )","(TC. |AG. )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TAG )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC |AAA )","(AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG |TGA )","(CG. )","(GG. )","(\S\S\S )");
        // Scenedesmus obliquus mitochondrial
        $triplets[22]=array("(TTT |TTC )","(TTA |TTG |CT. |TAG )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TCT |TCC |TCG |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TAA |TGA |TCA )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");
        // Thraustochytrium mitochondrial code
        $triplets[23]=array("(TTT |TTC )","(TTG |CT. )","(ATT |ATC |ATA )","(ATG )","(GT. )","(TC. |AGT |AGC )",
                        "(CC. )","(AC. )","(GC. )","(TAT |TAC )","(TTA |TAA |TAG |TGA )","(CAT |CAC )",
                        "(CAA |CAG )","(AAT |AAC )","(AAA |AAG )","(GAT |GAC )","(GAA |GAG )","(TGT |TGC )",
                        "(TGG )","(CG. |AGA |AGG )","(GG. )","(\S\S\S )");

        // place a space after each triplete in the sequence
        $temp = chunk_split($seq,3,' ');

        // replace triplets by corresponding amnoacid
        $peptide = preg_replace ($triplets[$genetic_code], $aminoacids, $temp);

        // return peptide sequence
        return trim($peptide);
}

// #################################################################################
// Section: translate_DNA_to_protein_customcode() function
// #################################################################################
function translate_DNA_to_protein_customcode($seq,$gc){
        // More info: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode

        // The sequence is chopped and @ is inserted after each triplete
        $temp=chunk_split($seq,3,' ');

        // each triplete replace by corresponding amnoacid
        $temp = str_replace ("TTT ",substr($gc, 0, 1)."  ",$temp);
        $temp = str_replace ("TTC ",substr($gc, 1, 1)."  ",$temp);
        $temp = str_replace ("TTA ",substr($gc, 2, 1)."  ",$temp);
        $temp = str_replace ("TTG ",substr($gc, 3, 1)."  ",$temp);
        $temp = str_replace ("TCT ",substr($gc, 4, 1)."  ",$temp);
        $temp = str_replace ("TCC ",substr($gc, 5, 1)."  ",$temp);
        $temp = str_replace ("TCA ",substr($gc, 6, 1)."  ",$temp);
        $temp = str_replace ("TCG ",substr($gc, 7, 1)."  ",$temp);
        $temp = str_replace ("TAT ",substr($gc, 8, 1)."  ",$temp);
        $temp = str_replace ("TAC ",substr($gc, 9, 1)."  ",$temp);
        $temp = str_replace ("TAA ",substr($gc,10, 1)."  ",$temp);
        $temp = str_replace ("TAG ",substr($gc,11, 1)."  ",$temp);
        $temp = str_replace ("TGT ",substr($gc,12, 1)."  ",$temp);
        $temp = str_replace ("TGC ",substr($gc,13, 1)."  ",$temp);
        $temp = str_replace ("TGA ",substr($gc,14, 1)."  ",$temp);
        $temp = str_replace ("TGG ",substr($gc,15, 1)."  ",$temp);
        $temp = str_replace ("CTT ",substr($gc,16, 1)."  ",$temp);
        $temp = str_replace ("CTC ",substr($gc,17, 1)."  ",$temp);
        $temp = str_replace ("CTA ",substr($gc,18, 1)."  ",$temp);
        $temp = str_replace ("CTG ",substr($gc,19, 1)."  ",$temp);
        $temp = str_replace ("CCT ",substr($gc,20, 1)."  ",$temp);
        $temp = str_replace ("CCC ",substr($gc,21, 1)."  ",$temp);
        $temp = str_replace ("CCA ",substr($gc,22, 1)."  ",$temp);
        $temp = str_replace ("CCG ",substr($gc,23, 1)."  ",$temp);
        $temp = str_replace ("CAT ",substr($gc,24, 1)."  ",$temp);
        $temp = str_replace ("CAC ",substr($gc,25, 1)."  ",$temp);
        $temp = str_replace ("CAA ",substr($gc,26, 1)."  ",$temp);
        $temp = str_replace ("CAG ",substr($gc,27, 1)."  ",$temp);
        $temp = str_replace ("CGT ",substr($gc,28, 1)."  ",$temp);
        $temp = str_replace ("CGC ",substr($gc,29, 1)."  ",$temp);
        $temp = str_replace ("CGA ",substr($gc,30, 1)."  ",$temp);
        $temp = str_replace ("CGG ",substr($gc,31, 1)."  ",$temp);
        $temp = str_replace ("ATT ",substr($gc,32, 1)."  ",$temp);
        $temp = str_replace ("ATC ",substr($gc,33, 1)."  ",$temp);
        $temp = str_replace ("ATA ",substr($gc,34, 1)."  ",$temp);
        $temp = str_replace ("ATG ",substr($gc,35, 1)."  ",$temp);
        $temp = str_replace ("ACT ",substr($gc,36, 1)."  ",$temp);
        $temp = str_replace ("ACC ",substr($gc,37, 1)."  ",$temp);
        $temp = str_replace ("ACA ",substr($gc,38, 1)."  ",$temp);
        $temp = str_replace ("ACG ",substr($gc,39, 1)."  ",$temp);
        $temp = str_replace ("AAT ",substr($gc,40, 1)."  ",$temp);
        $temp = str_replace ("AAC ",substr($gc,41, 1)."  ",$temp);
        $temp = str_replace ("AAA ",substr($gc,42, 1)."  ",$temp);
        $temp = str_replace ("AAG ",substr($gc,43, 1)."  ",$temp);
        $temp = str_replace ("AGT ",substr($gc,44, 1)."  ",$temp);
        $temp = str_replace ("AGC ",substr($gc,45, 1)."  ",$temp);
        $temp = str_replace ("AGA ",substr($gc,46, 1)."  ",$temp);
        $temp = str_replace ("AGG ",substr($gc,47, 1)."  ",$temp);
        $temp = str_replace ("GTT ",substr($gc,48, 1)."  ",$temp);
        $temp = str_replace ("GTC ",substr($gc,49, 1)."  ",$temp);
        $temp = str_replace ("GTA ",substr($gc,50, 1)."  ",$temp);
        $temp = str_replace ("GTG ",substr($gc,51, 1)."  ",$temp);
        $temp = str_replace ("GCT ",substr($gc,52, 1)."  ",$temp);
        $temp = str_replace ("GCC ",substr($gc,53, 1)."  ",$temp);
        $temp = str_replace ("GCA ",substr($gc,54, 1)."  ",$temp);
        $temp = str_replace ("GCG ",substr($gc,55, 1)."  ",$temp);
        $temp = str_replace ("GAT ",substr($gc,56, 1)."  ",$temp);
        $temp = str_replace ("GAC ",substr($gc,57, 1)."  ",$temp);
        $temp = str_replace ("GAA ",substr($gc,58, 1)."  ",$temp);
        $temp = str_replace ("GAG ",substr($gc,59, 1)."  ",$temp);
        $temp = str_replace ("GGT ",substr($gc,60, 1)."  ",$temp);
        $temp = str_replace ("GGC ",substr($gc,61, 1)."  ",$temp);
        $temp = str_replace ("GGA ",substr($gc,62, 1)."  ",$temp);
        $temp = str_replace ("GGG ",substr($gc,63, 1)."  ",$temp);
        // no matching triplets -> X
        $temp = preg_replace ("(\S\S\S )", "X  ", $temp);
        $temp = substr ($temp, 0, strlen($r)-2);
        $prot = preg_replace ("/ /","",$temp);
        return $prot;
}

// #################################################################################
// Section: str_is_int() function
// #################################################################################
function str_is_int($str) {
        $var=intval($str);
        return ("$str"=="$var");
}
// #################################################################################
// Section: str_is_int() function
//
// code modified from perl degen v1.4 script on www.phylotools.com with permission from
// creator A.Swick. 
// Users is asked to cite the following:
// Zwick, A., Regier, J.C. & Zwickl, D.J. (2012). "Resolving Discrepancy between Nucleotides and Amino Acids in Deep-Level Arthropod Phylogenomics: Differentiating Serine Codons in 21-Amino-Acid Models". PLoS ONE 7(11): e47450.
// Regier, J.C., Shultz, J.W., Zwick, A., Hussey, A., Ball, B., Wetzer, R. Martin, J.W. & Cunningham, C.W. (2010). "Arthropod relationships revealed by phylogenomic analysis of nuclear protein-coding sequences". Nature 463: 1079-1083.
// #################################################################################
function degen_coding($seq,$genetic_code,$rframe){
		//input: $seq is the sequence, $ is the genetic_code and $rframe is the readingframe.
		
        // $deg_pat is the array containning the pattern in genetic codes
		// $deg_code is the output degenerated triplet code
        // Info has been extracted from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode
// standard genetic code
$deg_code[1][]="YTN"; $deg_pat[1][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN )"; # Leu, Phe & Leu
$deg_code[1][]="TTY"; $deg_pat[1][]="(TTT |TTC |TTY )"; # Phe
$deg_code[1][]="MGN"; $deg_pat[1][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN )"; # Arg, Arg & Ser2
$deg_code[1][]="ATH"; $deg_pat[1][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW )"; # Ile
$deg_code[1][]="ATN"; $deg_pat[1][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN )"; # Ile & Met
$deg_code[1][]="ATG"; $deg_pat[1][]="(ATG  )"; # Met
$deg_code[1][]="GTN"; $deg_pat[1][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN )"; # Val
$deg_code[1][]="TCN"; $deg_pat[1][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN )"; # Ser1
$deg_code[1][]="AGY"; $deg_pat[1][]="(AGT |AGC |AGY )"; # Ser2
$deg_code[1][]="CCN"; $deg_pat[1][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN )"; # Pro
$deg_code[1][]="ACN"; $deg_pat[1][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN )"; # Thr
$deg_code[1][]="GCN"; $deg_pat[1][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN )"; # Ala
$deg_code[1][]="TAY"; $deg_pat[1][]="(TAT |TAC |TAY )"; # Tyr
$deg_code[1][]="CAY"; $deg_pat[1][]="(CAT |CAC |CAY )"; # His
$deg_code[1][]="CAN"; $deg_pat[1][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN )"; # His & Gln
$deg_code[1][]="CAR"; $deg_pat[1][]="(CAA |CAG |CAR )"; # Gln
$deg_code[1][]="AAY"; $deg_pat[1][]="(AAT |AAC |AAY )"; # Asn
$deg_code[1][]="AAN"; $deg_pat[1][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN )"; # Asn & Lys
$deg_code[1][]="AAR"; $deg_pat[1][]="(AAA |AAG |AAR )"; # Lys
$deg_code[1][]="GAY"; $deg_pat[1][]="(GAT |GAC |GAY )"; # Asp
$deg_code[1][]="GAN"; $deg_pat[1][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN )"; # Asp & Glu
$deg_code[1][]="GAR"; $deg_pat[1][]="(GAA |GAG |GAR )"; # Glu
$deg_code[1][]="TGY"; $deg_pat[1][]="(TGT |TGC |TGY )"; # Cys
$deg_code[1][]="TGB"; $deg_pat[1][]="(TGK |TGS |TGB )"; # Cys & Trp
$deg_code[1][]="TGG"; $deg_pat[1][]="(TGG )"; # Trp
$deg_code[1][]="GGN"; $deg_pat[1][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN )"; # Gly
$deg_code[1][]="NNN"; $deg_pat[1][]="(NNN )"; # degen
$deg_code[1][]='---'; $deg_pat[1][]="(--- )"; # indel
$deg_code[1][]='???'; $deg_pat[1][]="(\S\S\S )"; # empty

//degen2 |# 2. The |Vertebrate |Mitochondrial |Code
            #        Code |2          Standard
            # AGA |   Ter | *          Arg | R
            # AGG |   Ter | *          Arg | R
            # AUA |   Met | M          Ile | I
            # UGA |   Trp | W          Ter | *
$deg_code[2][]="YTN"; $deg_pat[2][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[2][]="TTY"; $deg_pat[2][]="(TTT |TTC |TTY)"; # Phe
$deg_code[2][]="CGN"; $deg_pat[2][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN)"; # Arg
$deg_code[2][]="AGY"; $deg_pat[2][]="(AGC |AGT |AGY)"; # Ser2
$deg_code[2][]="ATY"; $deg_pat[2][]="(ATT |ATC |ATY)"; # Ile
$deg_code[2][]="ATN"; $deg_pat[2][]="(ATM |ATW |ATS |ATK |ATV |ATH |ATD |ATB |ATN)"; # Ile |& Met
$deg_code[2][]="ATR"; $deg_pat[2][]="(ATG |ATA |ATR)"; # Met
$deg_code[2][]="GTN"; $deg_pat[2][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[2][]="TCN"; $deg_pat[2][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[2][]="CCN"; $deg_pat[2][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[2][]="ACN"; $deg_pat[2][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[2][]="GCN"; $deg_pat[2][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[2][]="TAY"; $deg_pat[2][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[2][]="CAY"; $deg_pat[2][]="(CAT |CAC |CAY)"; # His
$deg_code[2][]="CAN"; $deg_pat[2][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[2][]="CAR"; $deg_pat[2][]="(CAA |CAG |CAR)"; # Gln
$deg_code[2][]="AAY"; $deg_pat[2][]="(AAT |AAC |AAY)"; # Asn
$deg_code[2][]="AAN"; $deg_pat[2][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[2][]="AAR"; $deg_pat[2][]="(AAA |AAG |AAR)"; # Lys
$deg_code[2][]="GAY"; $deg_pat[2][]="(GAT |GAC |GAY)"; # Asp
$deg_code[2][]="GAN"; $deg_pat[2][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[2][]="GAR"; $deg_pat[2][]="(GAA |GAG |GAR)"; # Glu
$deg_code[2][]="TGY"; $deg_pat[2][]="(TGT |TGC |TGY)"; # Cys
$deg_code[2][]="TGN"; $deg_pat[2][]="(TGM |TGW |TGS |TGK |TGV |TGH |TGD |TGB |TGN)"; # Cys |& Trp
$deg_code[2][]="TGR"; $deg_pat[2][]="(TGG |TGA |TGR)"; # Trp
$deg_code[2][]="GGN"; $deg_pat[2][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[2][]="NNN"; $deg_pat[2][]="(NNN)"; # degen
$deg_code[2][]="---"; $deg_pat[2][]="(---)"; # indel
$deg_code[2][]='???'; $deg_pat[2][]="(\S\S\S )"; # empty


//degen3 |# 3. The |Yeast |Mitochondrial |Code
            #        Code |3          Standard
            # AUA |   Met | M          Ile | I
            # CUU |   Thr | T          Leu | L
            # CUC |   Thr | T          Leu | L
            # CUA |   Thr | T          Leu | L
            # CUG |   Thr | T          Leu | L
            # UGA |   Trp | W          Ter | *
            # CGA |   absent |         Arg | R
            # CGC |   absent |         Arg | R
$deg_code[3][]="TTR"; $deg_pat[3][]="(TTA |TTG |TTR)"; # Leu
$deg_code[3][]="TTN"; $deg_pat[3][]="(TTM |TTW |TTS |TTK |TTV |TTH |TTD |TTB |TTN)"; # Leu |& Phe
$deg_code[3][]="TTY"; $deg_pat[3][]="(TTT |TTC |TTY)"; # Phe
$deg_code[3][]="MGN"; $deg_pat[3][]="(CGT |CGG |CGR |CGY |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg, Arg |& Ser2
$deg_code[3][]="ATY"; $deg_pat[3][]="(ATT |ATC |ATY)"; # Ile
$deg_code[3][]="ATN"; $deg_pat[3][]="(ATM |ATW |ATS |ATK |ATV |ATH |ATD |ATB |ATN)"; # Ile |& Met
$deg_code[3][]="ATR"; $deg_pat[3][]="(ATG |ATA |ATR)"; # Met
$deg_code[3][]="GTN"; $deg_pat[3][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[3][]="TCN"; $deg_pat[3][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[3][]="AGY"; $deg_pat[3][]="(AGT |AGC |AGY)"; # Ser2
$deg_code[3][]="CCN"; $deg_pat[3][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[3][]="ACN"; $deg_pat[3][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr1
$deg_code[3][]="CTN"; $deg_pat[3][]="(CTT |CTA |CTC |CTG |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN)"; # Thr2
$deg_code[3][]="GCN"; $deg_pat[3][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[3][]="TAY"; $deg_pat[3][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[3][]="CAY"; $deg_pat[3][]="(CAT |CAC |CAY)"; # His
$deg_code[3][]="CAN"; $deg_pat[3][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[3][]="CAR"; $deg_pat[3][]="(CAA |CAG |CAR)"; # Gln
$deg_code[3][]="AAY"; $deg_pat[3][]="(AAT |AAC |AAY)"; # Asn
$deg_code[3][]="AAN"; $deg_pat[3][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[3][]="AAR"; $deg_pat[3][]="(AAA |AAG |AAR)"; # Lys
$deg_code[3][]="GAY"; $deg_pat[3][]="(GAT |GAC |GAY)"; # Asp
$deg_code[3][]="GAN"; $deg_pat[3][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[3][]="GAR"; $deg_pat[3][]="(GAA |GAG |GAR)"; # Glu
$deg_code[3][]="TGY"; $deg_pat[3][]="(TGT |TGC |TGY)"; # Cys
$deg_code[3][]="TGN"; $deg_pat[3][]="(TGM |TGW |TGS |TGK |TGV |TGH |TGD |TGB |TGN)"; # Cys |& Trp
$deg_code[3][]="TGR"; $deg_pat[3][]="(TGG |TGA |TGR)"; # Trp
$deg_code[3][]="GGN"; $deg_pat[3][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[3][]="NNN"; $deg_pat[3][]="(NNN)"; # degen
$deg_code[3][]="---"; $deg_pat[3][]="(---)"; # indel
$deg_code[3][]='???'; $deg_pat[3][]="(\S\S\S )"; # empty


//degen4 |# 4. The |Mold, Protozoan, and |Coelenterate |Mitochondrial |Code |and |the |Mycoplasma/Spiroplasma |Code
            #        Code |4         Standard
            # UGA |   Trp | W          Ter | *
$deg_code[4][]="YTN"; $deg_pat[4][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[4][]="TTY"; $deg_pat[4][]="(TTT |TTC |TTY)"; # Phe
$deg_code[4][]="MGN"; $deg_pat[4][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg, Arg |& Ser2
$deg_code[4][]="ATH"; $deg_pat[4][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code[4][]="ATN"; $deg_pat[4][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code[4][]="ATG"; $deg_pat[4][]="(ATG)"; # Met
$deg_code[4][]="GTN"; $deg_pat[4][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[4][]="TCN"; $deg_pat[4][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[4][]="AGY"; $deg_pat[4][]="(AGT |AGC |AGY)"; # Ser2
$deg_code[4][]="CCN"; $deg_pat[4][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[4][]="ACN"; $deg_pat[4][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[4][]="GCN"; $deg_pat[4][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[4][]="TAY"; $deg_pat[4][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[4][]="CAY"; $deg_pat[4][]="(CAT |CAC |CAY)"; # His
$deg_code[4][]="CAN"; $deg_pat[4][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[4][]="CAR"; $deg_pat[4][]="(CAA |CAG |CAR)"; # Gln
$deg_code[4][]="AAY"; $deg_pat[4][]="(AAT |AAC |AAY)"; # Asn
$deg_code[4][]="AAN"; $deg_pat[4][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[4][]="AAR"; $deg_pat[4][]="(AAA |AAG |AAR)"; # Lys
$deg_code[4][]="GAY"; $deg_pat[4][]="(GAT |GAC |GAY)"; # Asp
$deg_code[4][]="GAN"; $deg_pat[4][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[4][]="GAR"; $deg_pat[4][]="(GAA |GAG |GAR)"; # Glu
$deg_code[4][]="TGY"; $deg_pat[4][]="(TGT |TGC |TGY)"; # Cys
$deg_code[4][]="TGB"; $deg_pat[4][]="(TGM |TGK |TGS |TGW |TGH |TGB |TGV |TGD |TGN)"; # Cys |& Trp
$deg_code[4][]="TGR"; $deg_pat[4][]="(TGA |TGG |TGR)"; # Trp
$deg_code[4][]="GGN"; $deg_pat[4][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[4][]="NNN"; $deg_pat[4][]="(NNN)"; # degen
$deg_code[4][]="---"; $deg_pat[4][]="(---)"; # indel
$deg_code[4][]='???'; $deg_pat[4][]="(\S\S\S )"; # empty


//degen5 |# 5. The |Invertebrate |Mitochondrial |Code
            #        Code |5          Standard
            # AGA |   Ser | S          Arg | R
            # AGG |   Ser | S          Arg | R
            # AUA |   Met | M          Ile | I
            # UGA |   Trp | W          Ter | *
$deg_code[5][]="YTN"; $deg_pat[5][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[5][]="TTY"; $deg_pat[5][]="(TTT |TTC |TTY)"; # Phe
$deg_code[5][]="CGN"; $deg_pat[5][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN)"; # Arg
$deg_code[5][]="AGN"; $deg_pat[5][]="(AGT |AGA |AGC |AGG |AGR |AGY |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN)"; # Ser2
$deg_code[5][]="ATY"; $deg_pat[5][]="(ATT |ATC |ATY)"; # Ile
$deg_code[5][]="ATN"; $deg_pat[5][]="(ATM |ATW |ATS |ATK |ATV |ATH |ATD |ATB |ATN)"; # Ile |& Met
$deg_code[5][]="ATR"; $deg_pat[5][]="(ATG |ATA |ATR)"; # Met
$deg_code[5][]="GTN"; $deg_pat[5][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[5][]="TCN"; $deg_pat[5][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[5][]="CCN"; $deg_pat[5][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[5][]="ACN"; $deg_pat[5][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[5][]="GCN"; $deg_pat[5][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[5][]="TAY"; $deg_pat[5][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[5][]="CAY"; $deg_pat[5][]="(CAT |CAC |CAY)"; # His
$deg_code[5][]="CAN"; $deg_pat[5][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[5][]="CAR"; $deg_pat[5][]="(CAA |CAG |CAR)"; # Gln
$deg_code[5][]="AAY"; $deg_pat[5][]="(AAT |AAC |AAY)"; # Asn
$deg_code[5][]="AAN"; $deg_pat[5][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[5][]="AAR"; $deg_pat[5][]="(AAA |AAG |AAR)"; # Lys
$deg_code[5][]="GAY"; $deg_pat[5][]="(GAT |GAC |GAY)"; # Asp
$deg_code[5][]="GAN"; $deg_pat[5][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[5][]="GAR"; $deg_pat[5][]="(GAA |GAG |GAR)"; # Glu
$deg_code[5][]="TGY"; $deg_pat[5][]="(TGT |TGC |TGY)"; # Cys
$deg_code[5][]="TGN"; $deg_pat[5][]="(TGM |TGW |TGS |TGK |TGV |TGH |TGD |TGB |TGN)"; # Cys |& Trp
$deg_code[5][]="TGR"; $deg_pat[5][]="(TGG |TGA |TGR)"; # Trp
$deg_code[5][]="GGN"; $deg_pat[5][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[5][]="NNN"; $deg_pat[5][]="(NNN)"; # degen
$deg_code[5][]="---"; $deg_pat[5][]="(---)"; # indel
$deg_code[5][]='???'; $deg_pat[5][]="(\S\S\S )"; # empty


//degen6 |# 6. The |Ciliate, Dasycladacean |and |Hexamita |Nuclear |Code
            #          Code |6       Standard
            # UAA |     Gln | Q        Ter | *
            # UAG |     Gln | Q        Ter | *
$deg_code[6][]="YTN"; $deg_pat[6][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[6][]="TTY"; $deg_pat[6][]="(TTT |TTC |TTY)"; # Phe
$deg_code[6][]="MGN"; $deg_pat[6][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg, Arg |& Ser2
$deg_code[6][]="ATH"; $deg_pat[6][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code[6][]="ATN"; $deg_pat[6][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code[6][]="ATG"; $deg_pat[6][]="(ATG)"; # Met
$deg_code[6][]="GTN"; $deg_pat[6][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[6][]="TCN"; $deg_pat[6][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[6][]="AGY"; $deg_pat[6][]="(AGT |AGC |AGY)"; # Ser2
$deg_code[6][]="CCN"; $deg_pat[6][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[6][]="ACN"; $deg_pat[6][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[6][]="GCN"; $deg_pat[6][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[6][]="YAR"; $deg_pat[6][]="(TAA |TAG |TAR |CAA |CAG |CAR |YAA |YAG |YAR)"; # Gln
$deg_code[6][]="TAY"; $deg_pat[6][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[6][]="CAY"; $deg_pat[6][]="(CAT |CAC |CAY)"; # His
$deg_code[6][]="YAN"; $deg_pat[6][]="(TAM |TAK |TAS |TAW |TAH |TAB |TAV |TAD |TAN |CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN |YAM |YAK |YAS |YAW |YAH |YAB |YAV |YAD |YAN)"; # Gln, Tyr |& His
$deg_code[6][]="AAY"; $deg_pat[6][]="(AAT |AAC |AAY)"; # Asn
$deg_code[6][]="AAN"; $deg_pat[6][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[6][]="AAR"; $deg_pat[6][]="(AAA |AAG |AAR)"; # Lys
$deg_code[6][]="GAY"; $deg_pat[6][]="(GAT |GAC |GAY)"; # Asp
$deg_code[6][]="GAN"; $deg_pat[6][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[6][]="GAR"; $deg_pat[6][]="(GAA |GAG |GAR)"; # Glu
$deg_code[6][]="TGY"; $deg_pat[6][]="(TGT |TGC |TGY)"; # Cys
$deg_code[6][]="TGB"; $deg_pat[6][]="(TGK |TGS |TGB)"; # Cys |& Trp
$deg_code[6][]="TGG"; $deg_pat[6][]="(TGG)"; # Trp
$deg_code[6][]="GGN"; $deg_pat[6][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[6][]="NNN"; $deg_pat[6][]="(NNN)"; # degen
$deg_code[6][]="---"; $deg_pat[6][]="(---)"; # indel
$deg_code[6][]='???'; $deg_pat[6][]="(\S\S\S )"; # empty


//degen9 |# 9. The |Echinoderm |and |Flatworm |Mitochondrial |Code
            #        Code |9        Standard
            # AAA |     Asn | N        Lys |K
            # AGA |     Ser | S        Arg |R
            # AGG |     Ser | S        Arg |R
            # UGA |     Trp | W        Ter |*
$deg_code[9][]="YTN"; $deg_pat[9][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[9][]="TTY"; $deg_pat[9][]="(TTT |TTC |TTY)"; # Phe
$deg_code[9][]="CGN"; $deg_pat[9][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN)"; # Arg
$deg_code[9][]="AGN"; $deg_pat[9][]="(AGT |AGA |AGC |AGG |AGR |AGY |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN)"; # Ser2
$deg_code[9][]="ATH"; $deg_pat[9][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code[9][]="ATN"; $deg_pat[9][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code[9][]="ATG"; $deg_pat[9][]="(ATG)"; # Met
$deg_code[9][]="GTN"; $deg_pat[9][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[9][]="TCN"; $deg_pat[9][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[9][]="CCN"; $deg_pat[9][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[9][]="ACN"; $deg_pat[9][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[9][]="GCN"; $deg_pat[9][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[9][]="TAY"; $deg_pat[9][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[9][]="CAY"; $deg_pat[9][]="(CAT |CAC |CAY)"; # His
$deg_code[9][]="CAN"; $deg_pat[9][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[9][]="CAR"; $deg_pat[9][]="(CAA |CAG |CAR)"; # Gln
$deg_code[9][]="AAH"; $deg_pat[9][]="(AAA |AAT |AAC |AAY |AAM |AAW |AAH)"; # Asn
$deg_code[9][]="AAN"; $deg_pat[9][]="(AAR |AAK |AAS |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[9][]="AAG"; $deg_pat[9][]="(AAG)"; # Lys
$deg_code[9][]="GAY"; $deg_pat[9][]="(GAT |GAC |GAY)"; # Asp
$deg_code[9][]="GAN"; $deg_pat[9][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[9][]="GAR"; $deg_pat[9][]="(GAA |GAG |GAR)"; # Glu
$deg_code[9][]="TGY"; $deg_pat[9][]="(TGT |TGC |TGY)"; # Cys
$deg_code[9][]="TGN"; $deg_pat[9][]="(TGM |TGW |TGS |TGK |TGV |TGH |TGD |TGB |TGN)"; # Cys |& Trp
$deg_code[9][]="TGR"; $deg_pat[9][]="(TGG |TGA |TGR)"; # Trp
$deg_code[9][]="GGN"; $deg_pat[9][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[9][]="NNN"; $deg_pat[9][]="(NNN)"; # degen
$deg_code[9][]="---"; $deg_pat[9][]="(---)"; # indel
$deg_code[9][]='???'; $deg_pat[9][]="(\S\S\S )"; # empty


//degen10 |# 10. The |Euplotid |Nuclear |Code
            #          Code |10     Standard
            # UGA |     Cys | C        Ter | *
$deg_code[10][]="YTN"; $deg_pat[10][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[10][]="TTY"; $deg_pat[10][]="(TTT |TTC |TTY)"; # Phe
$deg_code[10][]="MGN"; $deg_pat[10][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg, Arg |& Ser2
$deg_code[10][]="ATH"; $deg_pat[10][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code[10][]="ATN"; $deg_pat[10][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code[10][]="ATG"; $deg_pat[10][]="(ATG)"; # Met
$deg_code[10][]="GTN"; $deg_pat[10][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[10][]="TCN"; $deg_pat[10][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[10][]="AGY"; $deg_pat[10][]="(AGT |AGC |AGY)"; # Ser2
$deg_code[10][]="CCN"; $deg_pat[10][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[10][]="ACN"; $deg_pat[10][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[10][]="GCN"; $deg_pat[10][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[10][]="TAY"; $deg_pat[10][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[10][]="CAY"; $deg_pat[10][]="(CAT |CAC |CAY)"; # His
$deg_code[10][]="CAN"; $deg_pat[10][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[10][]="CAR"; $deg_pat[10][]="(CAA |CAG |CAR)"; # Gln
$deg_code[10][]="AAY"; $deg_pat[10][]="(AAT |AAC |AAY)"; # Asn
$deg_code[10][]="AAN"; $deg_pat[10][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[10][]="AAR"; $deg_pat[10][]="(AAA |AAG |AAR)"; # Lys
$deg_code[10][]="GAY"; $deg_pat[10][]="(GAT |GAC |GAY)"; # Asp
$deg_code[10][]="GAN"; $deg_pat[10][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[10][]="GAR"; $deg_pat[10][]="(GAA |GAG |GAR)"; # Glu
$deg_code[10][]="TGH"; $deg_pat[10][]="(TGT |TGC |TGA |TGH |TGY |TGM |TGW)"; # Cys
$deg_code[10][]="TGN"; $deg_pat[10][]="(TGR |TGK |TGS |TGB |TGV |TGD |TGN)"; # Cys |& Trp
$deg_code[10][]="TGG"; $deg_pat[10][]="(TGG)"; # Trp
$deg_code[10][]="GGN"; $deg_pat[10][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[10][]="NNN"; $deg_pat[10][]="(NNN)"; # degen
$deg_code[10][]="---"; $deg_pat[10][]="(---)"; # indel
$deg_code[10][]='???'; $deg_pat[10][]="(\S\S\S )"; # empty


//degen11 |# 11. The |Bacterial, Archaeal |and |Plant |Plastid |Code
$deg_code[11][]="YTN"; $deg_pat[11][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[11][]="TTY"; $deg_pat[11][]="(TTT |TTC |TTY)"; # Phe
$deg_code[11][]="MGN"; $deg_pat[11][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg, Arg |& Ser2
$deg_code[11][]="ATH"; $deg_pat[11][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code[11][]="ATN"; $deg_pat[11][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code[11][]="ATG"; $deg_pat[11][]="(ATG)"; # Met
$deg_code[11][]="GTN"; $deg_pat[11][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[11][]="TCN"; $deg_pat[11][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[11][]="AGY"; $deg_pat[11][]="(AGT |AGC |AGY)"; # Ser2
$deg_code[11][]="CCN"; $deg_pat[11][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[11][]="ACN"; $deg_pat[11][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[11][]="GCN"; $deg_pat[11][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[11][]="TAY"; $deg_pat[11][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[11][]="CAY"; $deg_pat[11][]="(CAT |CAC |CAY)"; # His
$deg_code[11][]="CAN"; $deg_pat[11][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[11][]="CAR"; $deg_pat[11][]="(CAA |CAG |CAR)"; # Gln
$deg_code[11][]="AAY"; $deg_pat[11][]="(AAT |AAC |AAY)"; # Asn
$deg_code[11][]="AAN"; $deg_pat[11][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[11][]="AAR"; $deg_pat[11][]="(AAA |AAG |AAR)"; # Lys
$deg_code[11][]="GAY"; $deg_pat[11][]="(GAT |GAC |GAY)"; # Asp
$deg_code[11][]="GAN"; $deg_pat[11][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[11][]="GAR"; $deg_pat[11][]="(GAA |GAG |GAR)"; # Glu
$deg_code[11][]="TGY"; $deg_pat[11][]="(TGT |TGC |TGY)"; # Cys
$deg_code[11][]="TGB"; $deg_pat[11][]="(TGK |TGS |TGB)"; # Cys |& Trp
$deg_code[11][]="TGG"; $deg_pat[11][]="(TGG)"; # Trp
$deg_code[11][]="GGN"; $deg_pat[11][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[11][]="NNN"; $deg_pat[11][]="(NNN)"; # degen
$deg_code[11][]="---"; $deg_pat[11][]="(---)"; # indel
$deg_code[11][]='???'; $deg_pat[11][]="(\S\S\S )"; # empty


//degen12 |# 12. The |Alternative |Yeast |Nuclear |Code
            #           Code |12      Standard
            # CUG |      Ser |         Leu
$deg_code[12][]="YTN"; $deg_pat[12][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu; Phe |& Leu; Phe, Leu |& Ser3
$deg_code[12][]="CTG"; $deg_pat[12][]="(CTG)"; # Ser3
$deg_code[12][]="TTY"; $deg_pat[12][]="(TTT |TTC |TTY)"; # Phe
$deg_code[12][]="MGN"; $deg_pat[12][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg; Arg |& Ser2
$deg_code[12][]="ATH"; $deg_pat[12][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code[12][]="ATN"; $deg_pat[12][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code[12][]="ATG"; $deg_pat[12][]="(ATG)"; # Met
$deg_code[12][]="GTN"; $deg_pat[12][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[12][]="TCN"; $deg_pat[12][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[12][]="AGY"; $deg_pat[12][]="(AGT |AGC |AGY)"; # Ser2
$deg_code[12][]="CCN"; $deg_pat[12][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[12][]="ACN"; $deg_pat[12][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[12][]="GCN"; $deg_pat[12][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[12][]="TAY"; $deg_pat[12][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[12][]="CAY"; $deg_pat[12][]="(CAT |CAC |CAY)"; # His
$deg_code[12][]="CAN"; $deg_pat[12][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[12][]="CAR"; $deg_pat[12][]="(CAA |CAG |CAR)"; # Gln
$deg_code[12][]="AAY"; $deg_pat[12][]="(AAT |AAC |AAY)"; # Asn
$deg_code[12][]="AAN"; $deg_pat[12][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[12][]="AAR"; $deg_pat[12][]="(AAA |AAG |AAR)"; # Lys
$deg_code[12][]="GAY"; $deg_pat[12][]="(GAT |GAC |GAY)"; # Asp
$deg_code[12][]="GAN"; $deg_pat[12][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[12][]="GAR"; $deg_pat[12][]="(GAA |GAG |GAR)"; # Glu
$deg_code[12][]="TGY"; $deg_pat[12][]="(TGT |TGC |TGY)"; # Cys
$deg_code[12][]="TGB"; $deg_pat[12][]="(TGK |TGS |TGB)"; # Cys |& Trp
$deg_code[12][]="TGG"; $deg_pat[12][]="(TGG)"; # Trp
$deg_code[12][]="GGN"; $deg_pat[12][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[12][]="NNN"; $deg_pat[12][]="(NNN)"; # degen
$deg_code[12][]="---"; $deg_pat[12][]="(---)"; # indel
$deg_code[12][]='???'; $deg_pat[12][]="(\S\S\S )"; # empty


//degen13 |# 13. The |Ascidian |Mitochondrial |Code
            #        Code |13         Standard
            # AGA |   Gly | G          Arg | R
            # AGG |   Gly | G          Arg | R
            # AUA |   Met | M          Ile | I
            # UGA |   Trp | W          Ter | *
$deg_code[13][]="YTN"; $deg_pat[13][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[13][]="TTY"; $deg_pat[13][]="(TTT |TTC |TTY)"; # Phe
$deg_code[13][]="CGN"; $deg_pat[13][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN)"; # Arg
$deg_code[13][]="AGY"; $deg_pat[13][]="(AGC |AGT |AGY)"; # Ser2
$deg_code[13][]="ATY"; $deg_pat[13][]="(ATT |ATC |ATY)"; # Ile
$deg_code[13][]="ATN"; $deg_pat[13][]="(ATM |ATW |ATS |ATK |ATV |ATH |ATD |ATB |ATN)"; # Ile |& Met
$deg_code[13][]="ATR"; $deg_pat[13][]="(ATG |ATA |ATR)"; # Met
$deg_code[13][]="GTN"; $deg_pat[13][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[13][]="TCN"; $deg_pat[13][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[13][]="CCN"; $deg_pat[13][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[13][]="ACN"; $deg_pat[13][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[13][]="GCN"; $deg_pat[13][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[13][]="TAY"; $deg_pat[13][]="(TAT |TAC |TAY)"; # Tyr
$deg_code[13][]="CAY"; $deg_pat[13][]="(CAT |CAC |CAY)"; # His
$deg_code[13][]="CAN"; $deg_pat[13][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[13][]="CAR"; $deg_pat[13][]="(CAA |CAG |CAR)"; # Gln
$deg_code[13][]="AAY"; $deg_pat[13][]="(AAT |AAC |AAY)"; # Asn
$deg_code[13][]="AAN"; $deg_pat[13][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[13][]="AAR"; $deg_pat[13][]="(AAA |AAG |AAR)"; # Lys
$deg_code[13][]="GAY"; $deg_pat[13][]="(GAT |GAC |GAY)"; # Asp
$deg_code[13][]="GAN"; $deg_pat[13][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[13][]="GAR"; $deg_pat[13][]="(GAA |GAG |GAR)"; # Glu
$deg_code[13][]="TGY"; $deg_pat[13][]="(TGT |TGC |TGY)"; # Cys
$deg_code[13][]="TGN"; $deg_pat[13][]="(TGM |TGW |TGS |TGK |TGV |TGH |TGD |TGB |TGN)"; # Cys |& Trp
$deg_code[13][]="TGR"; $deg_pat[13][]="(TGG |TGA |TGR)"; # Trp
$deg_code[13][]="RGN"; $deg_pat[13][]="(AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |AGA |AGG |AGR |GGA |GGC |GGG |GGT |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN |RGC |RGA |RGT |RGG |RGR |RGY |RGM |RGK |RGS |RGW |RGH |RGB |RGV |RGD |RGN)"; # Gly |& Ser2
$deg_code[13][]="NNN"; $deg_pat[13][]="(NNN)"; # degen
$deg_code[13][]="---"; $deg_pat[13][]="(---)"; # indel
$deg_code[13][]='???'; $deg_pat[13][]="(\S\S\S )"; # empty


//degen14 |# 14. The |Alternative |Flatworm |Mitochondrial |Code
            #          Code |14      Standard
            # AAA |     Asn | N       Lys | K
            # AGA |     Ser | S       Arg | R
            # AGG |     Ser | S       Arg | R
            # UAA |     Tyr | Y       Ter | *
            # UGA |     Trp | W       Ter | *
$deg_code[14][]="YTN"; $deg_pat[14][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code[14][]="TTY"; $deg_pat[14][]="(TTT |TTC |TTY)"; # Phe
$deg_code[14][]="CGN"; $deg_pat[14][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN)"; # Arg
$deg_code[14][]="AGN"; $deg_pat[14][]="(AGT |AGA |AGC |AGG |AGR |AGY |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN)"; # Ser2
$deg_code[14][]="ATH"; $deg_pat[14][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code[14][]="ATN"; $deg_pat[14][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code[14][]="ATG"; $deg_pat[14][]="(ATG)"; # Met
$deg_code[14][]="GTN"; $deg_pat[14][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code[14][]="TCN"; $deg_pat[14][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN)"; # Ser1
$deg_code[14][]="CCN"; $deg_pat[14][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code[14][]="ACN"; $deg_pat[14][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code[14][]="GCN"; $deg_pat[14][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code[14][]="TAH"; $deg_pat[14][]="(TAT |TAC |TAA |TAH |TAY |TAM |TAW)"; # Tyr
$deg_code[14][]="CAY"; $deg_pat[14][]="(CAT |CAC |CAY)"; # His
$deg_code[14][]="CAN"; $deg_pat[14][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code[14][]="CAR"; $deg_pat[14][]="(CAA |CAG |CAR)"; # Gln
$deg_code[14][]="AAH"; $deg_pat[14][]="(AAA |AAT |AAC |AAY |AAM |AAW |AAH)"; # Asn
$deg_code[14][]="AAN"; $deg_pat[14][]="(AAR |AAK |AAS |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code[14][]="AAG"; $deg_pat[14][]="(AAG)"; # Lys
$deg_code[14][]="GAY"; $deg_pat[14][]="(GAT |GAC |GAY)"; # Asp
$deg_code[14][]="GAN"; $deg_pat[14][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code[14][]="GAR"; $deg_pat[14][]="(GAA |GAG |GAR)"; # Glu
$deg_code[14][]="TGY"; $deg_pat[14][]="(TGT |TGC |TGY)"; # Cys
$deg_code[14][]="TGN"; $deg_pat[14][]="(TGM |TGW |TGS |TGK |TGV |TGH |TGD |TGB |TGN)"; # Cys |& Trp
$deg_code[14][]="TGR"; $deg_pat[14][]="(TGG |TGA |TGR)"; # Trp
$deg_code[14][]="GGN"; $deg_pat[14][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code[14][]="NNN"; $deg_pat[14][]="(NNN)"; # degen
$deg_code[14][]="---"; $deg_pat[14][]="(---)"; # indel
$deg_code[14][]='???'; $deg_pat[14][]="(\S\S\S )"; # empty


//degenS - standard code with Ser1 and Ser2 to Ser2
$deg_code['S'][]="YTN"; $deg_pat['S'][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code['S'][]="TTY"; $deg_pat['S'][]="(TTT |TTC |TTY)"; # Phe
$deg_code['S'][]="MGN"; $deg_pat['S'][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg, Arg |& Ser2
$deg_code['S'][]="ATH"; $deg_pat['S'][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code['S'][]="ATN"; $deg_pat['S'][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code['S'][]="ATG"; $deg_pat['S'][]="(ATG)"; # Met
$deg_code['S'][]="GTN"; $deg_pat['S'][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code['S'][]="AGY"; $deg_pat['S'][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN |AGT |AGC |AGY)"; # Ser1 |and |Ser2 |to Ser2
$deg_code['S'][]="CCN"; $deg_pat['S'][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code['S'][]="ACN"; $deg_pat['S'][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code['S'][]="GCN"; $deg_pat['S'][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code['S'][]="TAY"; $deg_pat['S'][]="(TAT |TAC |TAY)"; # Tyr
$deg_code['S'][]="CAY"; $deg_pat['S'][]="(CAT |CAC |CAY)"; # His
$deg_code['S'][]="CAN"; $deg_pat['S'][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code['S'][]="CAR"; $deg_pat['S'][]="(CAA |CAG |CAR)"; # Gln
$deg_code['S'][]="AAY"; $deg_pat['S'][]="(AAT |AAC |AAY)"; # Asn
$deg_code['S'][]="AAN"; $deg_pat['S'][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code['S'][]="AAR"; $deg_pat['S'][]="(AAA |AAG |AAR)"; # Lys
$deg_code['S'][]="GAY"; $deg_pat['S'][]="(GAT |GAC |GAY)"; # Asp
$deg_code['S'][]="GAN"; $deg_pat['S'][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code['S'][]="GAR"; $deg_pat['S'][]="(GAA |GAG |GAR)"; # Glu
$deg_code['S'][]="TGY"; $deg_pat['S'][]="(TGT |TGC |TGY)"; # Cys
$deg_code['S'][]="TGB"; $deg_pat['S'][]="(TGK |TGS |TGB)"; # Cys |& Trp
$deg_code['S'][]="TGG"; $deg_pat['S'][]="(TGG)"; # Trp
$deg_code['S'][]="GGN"; $deg_pat['S'][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code['S'][]="NNN"; $deg_pat['S'][]="(NNN)"; # degen
$deg_code['S'][]="---"; $deg_pat['S'][]="(---)"; # indel
$deg_code['S'][]='???'; $deg_pat['S'][]="(\S\S\S )"; # empty


//degenSZ - standard code with Ser1 and Ser2 to NNN
$deg_code['SZ'][]="YTN"; $deg_pat['SZ'][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code['SZ'][]="TTY"; $deg_pat['SZ'][]="(TTT |TTC |TTY)"; # Phe
$deg_code['SZ'][]="MGN"; $deg_pat['SZ'][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg, Arg |& Ser2
$deg_code['SZ'][]="ATH"; $deg_pat['SZ'][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code['SZ'][]="ATN"; $deg_pat['SZ'][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code['SZ'][]="ATG"; $deg_pat['SZ'][]="(ATG)"; # Met
$deg_code['SZ'][]="GTN"; $deg_pat['SZ'][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code['SZ'][]="CCN"; $deg_pat['SZ'][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code['SZ'][]="ACN"; $deg_pat['SZ'][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code['SZ'][]="GCN"; $deg_pat['SZ'][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code['SZ'][]="TAY"; $deg_pat['SZ'][]="(TAT |TAC |TAY)"; # Tyr
$deg_code['SZ'][]="CAY"; $deg_pat['SZ'][]="(CAT |CAC |CAY)"; # His
$deg_code['SZ'][]="CAN"; $deg_pat['SZ'][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code['SZ'][]="CAR"; $deg_pat['SZ'][]="(CAA |CAG |CAR)"; # Gln
$deg_code['SZ'][]="AAY"; $deg_pat['SZ'][]="(AAT |AAC |AAY)"; # Asn
$deg_code['SZ'][]="AAN"; $deg_pat['SZ'][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code['SZ'][]="AAR"; $deg_pat['SZ'][]="(AAA |AAG |AAR)"; # Lys
$deg_code['SZ'][]="GAY"; $deg_pat['SZ'][]="(GAT |GAC |GAY)"; # Asp
$deg_code['SZ'][]="GAN"; $deg_pat['SZ'][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code['SZ'][]="GAR"; $deg_pat['SZ'][]="(GAA |GAG |GAR)"; # Glu
$deg_code['SZ'][]="TGY"; $deg_pat['SZ'][]="(TGT |TGC |TGY)"; # Cys
$deg_code['SZ'][]="TGB"; $deg_pat['SZ'][]="(TGK |TGS |TGB)"; # Cys |& Trp
$deg_code['SZ'][]="TGG"; $deg_pat['SZ'][]="(TGG)"; # Trp
$deg_code['SZ'][]="GGN"; $deg_pat['SZ'][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code['SZ'][]="NNN"; $deg_pat['SZ'][]="(NNN |TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN |AGT |AGC |AGY)"; # degen
$deg_code['SZ'][]="---"; $deg_pat['SZ'][]="(---)"; # indel
$deg_code['SZ'][]='???'; $deg_pat['SZ'][]="(\S\S\S )"; # empty


//degenZ - standard code with Ser1 and Ser2 to Ser1
$deg_code['Z'][]="YTN"; $deg_pat['Z'][]="(TTM |TTK |TTS |TTW |TTH |TTB |TTV |TTD |TTN |TTA |TTG |TTR |CTA |CTC |CTG |CTT |CTR |CTY |CTM |CTK |CTS |CTW |CTH |CTB |CTV |CTD |CTN |YTC |YTA |YTT |YTG |YTR |YTY |YTM |YTK |YTS |YTW |YTH |YTB |YTV |YTD |YTN)"; # Leu, Phe |& Leu
$deg_code['Z'][]="TTY"; $deg_pat['Z'][]="(TTT |TTC |TTY)"; # Phe
$deg_code['Z'][]="MGN"; $deg_pat['Z'][]="(CGT |CGA |CGC |CGG |CGR |CGY |CGM |CGK |CGS |CGW |CGH |CGB |CGV |CGD |CGN |AGA |AGG |AGR |AGM |AGK |AGS |AGW |AGH |AGB |AGV |AGD |AGN |MGC |MGA |MGT |MGG |MGR |MGY |MGM |MGK |MGS |MGW |MGH |MGB |MGV |MGD |MGN)"; # Arg, Arg |& Ser2
$deg_code['Z'][]="ATH"; $deg_pat['Z'][]="(ATT |ATC |ATA |ATH |ATY |ATM |ATW)"; # Ile
$deg_code['Z'][]="ATN"; $deg_pat['Z'][]="(ATR |ATK |ATS |ATB |ATV |ATD |ATN)"; # Ile |& Met
$deg_code['Z'][]="ATG"; $deg_pat['Z'][]="(ATG)"; # Met
$deg_code['Z'][]="GTN"; $deg_pat['Z'][]="(GTT |GTA |GTC |GTG |GTR |GTY |GTM |GTK |GTS |GTW |GTH |GTB |GTV |GTD |GTN)"; # Val
$deg_code['Z'][]="TCN"; $deg_pat['Z'][]="(TCT |TCA |TCC |TCG |TCR |TCY |TCM |TCK |TCS |TCW |TCH |TCB |TCV |TCD |TCN |AGT |AGC |AGY)"; # Ser1 |and |Ser2 |to Ser1
$deg_code['Z'][]="CCN"; $deg_pat['Z'][]="(CCT |CCA |CCC |CCG |CCR |CCY |CCM |CCK |CCS |CCW |CCH |CCB |CCV |CCD |CCN)"; # Pro
$deg_code['Z'][]="ACN"; $deg_pat['Z'][]="(ACT |ACA |ACC |ACG |ACR |ACY |ACM |ACK |ACS |ACW |ACH |ACB |ACV |ACD |ACN)"; # Thr
$deg_code['Z'][]="GCN"; $deg_pat['Z'][]="(GCT |GCA |GCC |GCG |GCR |GCY |GCM |GCK |GCS |GCW |GCH |GCB |GCV |GCD |GCN)"; # Ala
$deg_code['Z'][]="TAY"; $deg_pat['Z'][]="(TAT |TAC |TAY)"; # Tyr
$deg_code['Z'][]="CAY"; $deg_pat['Z'][]="(CAT |CAC |CAY)"; # His
$deg_code['Z'][]="CAN"; $deg_pat['Z'][]="(CAM |CAK |CAS |CAW |CAH |CAB |CAV |CAD |CAN)"; # His |& Gln
$deg_code['Z'][]="CAR"; $deg_pat['Z'][]="(CAA |CAG |CAR)"; # Gln
$deg_code['Z'][]="AAY"; $deg_pat['Z'][]="(AAT |AAC |AAY)"; # Asn
$deg_code['Z'][]="AAN"; $deg_pat['Z'][]="(AAM |AAK |AAS |AAW |AAH |AAB |AAV |AAD |AAN)"; # Asn |& Lys
$deg_code['Z'][]="AAR"; $deg_pat['Z'][]="(AAA |AAG |AAR)"; # Lys
$deg_code['Z'][]="GAY"; $deg_pat['Z'][]="(GAT |GAC |GAY)"; # Asp
$deg_code['Z'][]="GAN"; $deg_pat['Z'][]="(GAM |GAK |GAS |GAW |GAH |GAB |GAV |GAD |GAN)"; # Asp |& Glu
$deg_code['Z'][]="GAR"; $deg_pat['Z'][]="(GAA |GAG |GAR)"; # Glu
$deg_code['Z'][]="TGY"; $deg_pat['Z'][]="(TGT |TGC |TGY)"; # Cys
$deg_code['Z'][]="TGB"; $deg_pat['Z'][]="(TGK |TGS |TGB)"; # Cys |& Trp
$deg_code['Z'][]="TGG"; $deg_pat['Z'][]="(TGG)"; # Trp
$deg_code['Z'][]="GGN"; $deg_pat['Z'][]="(GGT |GGA |GGC |GGG |GGR |GGY |GGM |GGK |GGS |GGW |GGH |GGB |GGV |GGD |GGN)"; # Gly
$deg_code['Z'][]="NNN"; $deg_pat['Z'][]="(NNN)"; # degen
$deg_code['Z'][]="---"; $deg_pat['Z'][]="(---)"; # indel
$deg_code['Z'][]='???'; $deg_pat['Z'][]="(\S\S\S )"; # empty



		$degen_seq = '';
		$deg_seq_end = '';
		//fix trimmings
		if ($rframe == "2"){ $seqstr = substr($seq, 1);$degen_seq .= substr($seq, 0,1);}
		if ($rframe == "3"){ $seqstr = substr($seq, 2);$degen_seq .= substr($seq, 0,2);}
		
		// place a space after each triplete in the sequence
		$temp = chunk_split($seqstr,3,' ');
		
		if (substr($temp, -1, 1) != " "){
			for ($i=1;$i<=3;$i++){
				if (substr($temp, -1, 1) != " "){$deg_seq_end = substr($temp, -1, 1) + $deg_seq_end; $temp = substr($temp, 0, -1);}
			}
		}

		// replace triplets by corresponding degen code
		$seq_rep = preg_replace ($deg_pat[$genetic_code], $deg_code[$genetic_code], $temp);
		
		//
		//paste together resulting seq
		$degen_seq = $degen_seq . $seq_rep . $deg_seq_end;
        // return peptide sequence
        return trim($degen_seq);
}
?>
