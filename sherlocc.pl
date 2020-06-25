#!/usr/bin/perl -w

%hydro = ("I", 1, "V", 0.967, "L", 0.922, "F", 0.811, "C", 0.778, "M", 0.711, "A", 0.700, "G", 0.456, "T", 0.422, "W", 0.400, "S", 0.411, "Y", 0.356, "P", 0.322, "H", 0.144, "E", 0.111, "Q", 0.111, "D", 0.111, "N", 0.111, "K", 0.067, "R", 0);

%charge = ("I", 0, "V", 0, "L", 0, "F", 0, "C", 0, "M", 0, "A", 0, "G", 0, "T", 0, "W", 0, "S", 0, "Y", 0, "P", 0, "H", 0, "E", -1, "Q", 0, "D", -1, "N", 0, "K", 1, "R", 1);

#By Matthieu Chartier

#Ask to define foldindex window size
print "\nStart at stage 3? y or n\nStage 2 files are needed and are available at bcb.med.usherbrooke.ca/pages/sherlocc\nIf y, pre-existing files in directory stage2/ will be deleted\n";
chomp($start_at_s3 = <STDIN>);
$start_at_s3="n" if($start_at_s3 eq "");

if($start_at_s3 eq "n"){
    cleanup("stage1") if(-d "stage1");
    cleanup("stage2") if(-d "stage2");
    mkdir "stage1";
    mkdir "stage2";
    print "\nMake sure there are only extracted Pfam alignment files in pfam_files/\n";
}elsif($start_at_s3=~/^y$/i){
    cleanup("stage2") if(-d "stage2");
    mkdir "stage2";
    print "\n! Please copy stage2 files in the newly created stage2/ directory!\n";
}

#Define rareness threshold
print "\nEnter Sherlocc threshold:\n";
chomp($avrg_lim = <STDIN>);

#Ask to define foldindex window size
print "\nEnter a window size for unstructured region prediction (default: 51):\n";
chomp($win_size = <STDIN>);
if($win_size eq ""){
    $win_size=51;
}

#Ask to cleanup stage 1 and stage 2 files
if($start_at_s3 eq "n"){
    print "\nDo you want to clean up stage 1 and stage 2 files after the analysis (y or n):\n";
    chomp($cleanup_dirs = <STDIN>);
    $cleanup_dirs="n" if($cleanup_dirs eq "");
}else{
    $cleanup_dirs="n";
}


#Ask if user wants to create HTML outputs
print "\nCreate html outputs? May take disk space if many pfam files (y/n):\n";
chomp($cr_html = <STDIN>);
unless($cr_html =~ /^y$/i){
    $cr_html="n";
}

#Ask for a dir name to store results
print "\nEnter a directory name in which to store results (in results/):\n";
chomp($store_res_dir = <STDIN>);

while(){
    if(-d "results/$store_res_dir"){
	print "\nDir $store_res_dir already exists.\nEmpty its content (a) or give a new name(b)?\n";
	chomp($del_results_dir = <STDIN>);
	if($del_results_dir eq "a"){
	    cleanup("results/$store_res_dir");
	    mkdir "results/$store_res_dir" or die "cant create dir results/$store_res_dir";
	    last;
	}
	if($del_results_dir eq "b"){
	    print "\nEnter a new directory name to store results:\n";
	    chomp($store_res_dir = <STDIN>);
	}
    }else{
	mkdir "results/$store_res_dir" or die "cant create dir results/$store_res_dir";
	last;
    }
}

#Build translation table hash
open TRANSL_TBL, "<files/transl_tbl" or die "Cant open transl_tbl file";
foreach(<TRANSL_TBL>){
    @elements=split(/,/, $_);
}
foreach(@elements){
    $_=~s/^\s+//;
    $transl_tbl{$1}=$2 if($_=~/'([a-z0-9\*]+)'=\>'([a-z0-9\*]+)'/i);
}
print "\nDone building translation table\n";

#Open codon usage file and store the table as an array
open CUT, "files/cut" or die "\nCan't open file codon usage table file\n";
foreach(<CUT>){
    push @cut, $_;
}
close CUT;

unless($start_at_s3 eq "y"){
    print "\n--------------- STAGE 1 ---------------\n";
    print "Retrieving info for each Pfam\n";
    @pfam_auto_files_list = glob "pfam_files/*";
    foreach(@pfam_auto_files_list){
	chomp($pfam_file = $_);
	$pfam_file =~ s/pfam_files\///;
	print "\n$pfam_file\n";
	
	open CONTENT, "pfam_files/$pfam_file" or die "\nCan't open PFAM alignment file $pfam_file\n";
	foreach(<CONTENT>){
	    push @pfam_file, $_;
	}
	$nb_index = $#pfam_file;
	
	#Open file in which to write data
	open ALIGN, ">stage1/pfam_$pfam_file" or die "Cant open output file in pfam_output";
	
	#Search the pfam alignment for the line at which alignment starts
	for($i=0; $i<=$nb_index; $i++){
	    if($pfam_file[$i] =~ /^[^#]/){
		$index = $i;
		last;
	    }
	}
	
	###########################################################
	#-----| FOR EACH PROTEIN IN THE PFAM ALIGNMENT FILE |-----#
	#-----|-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-|-----#
	###########################################################
	
	$protein_count=0;
	$protein_stored=0;
	$protein_no_cu=0;
	$protein_no_ebi=0;
	$protein_aa_not_matched=0;
	$protein_wrong_codon_nb=0;
	$protein_nuc_no_match=0;
	
      PROTEIN: for($z=$index; $z<=$nb_index; $z++){
	  if($pfam_file[$z] =~ /^([A-Z0-9]+_?[0-9A-Za-z]*)\//){
	      $uniprot_id = $1;
	      $protein_count++;
	      
	      #Get content from the uniprot website for this protein
	      $uniprot_url = "http://www.uniprot.org/uniprot/$uniprot_id";
	      use LWP::Simple;
	      @uniprot_content = ();
	      @split_uniprot_content = ();
	      @uniprot_content = get $uniprot_url or die "\nCan't get content from $uniprot_url";	
	      @split_uniprot_content = split(/\n/, "@uniprot_content");
	      
	      #Get organism's taxonomic ID	
	      foreach(@split_uniprot_content) {
		  if($_ =~ /http:\/\/www.ncbi.nlm.nih.gov\/Taxonomy\/Browser\/wwwtax.cgi\?lvl=0&amp;id=(\d+)"/i) {
		      $org_id = $1;
		      last;
		  }
	      }
	      
	      #Get Uniprot identifier
	      foreach(@split_uniprot_content){
		  $uniprot_id = $1 if($_ =~ /\<strong\>\s+([0-9a-z]+)\<\/strong\>\s+\($uniprot_id\)/i);
	      }
	      
	      $matched="0";
	      #Check if codon usage table of this organism is available
	      foreach(@cut){
		  if(/=($org_id)\|/){
		      $matched="1";
		      #print "\nCodon usage found";
		      last;
		  }
	      }
	      if($matched=="0"){ #Error #1
		  #print "\nCodon usage not available\n";
		  $protein_no_cu++;
		  push @protein_no_cu, $uniprot_id;
		  next PROTEIN;
	      }
	      
	      if($matched == "1"){
		  #Get pfam_ori and pfam_aa
		  chomp($pfam_align=$pfam_file[$z]);
		  $position=rindex($pfam_align, " ");
		  $pfam_aa_ori=substr $pfam_align, ($position + 1);
		  $pfam_aa_ori=~tr/a-z/A-Z/;
		  $pfam_aa=substr $pfam_align, ($position + 1);
		  $pfam_aa=~s/\.//g;
		  $pfam_aa=~s/-//g;
		  $pfam_aa=~tr/a-z/A-Z/;
		  
		  #Find ebi url from uniprot website
		  $ebi_url="";
		  foreach(@split_uniprot_content){
		      if($_=~/\<a\s+class="embl_cds"\s+href="(http:\/\/www\.ebi\.ac\.uk\/ena\/data\/view\/[\w]+)"\s+/i){
			  $ebi_url=$1;
			  last;
		      }
		  }
		  
		  if($ebi_url eq ""){ #Error #2
		      #print "\nCouldn't find cds from ebi.";
		      $protein_no_ebi++;
		      push @protein_no_ebi, $uniprot_id;
		      next PROTEIN;
		  }
		  $ebi_url.="&display=text";
		  #print "\n$ebi_url";
		  
		  ############################################################
		  #########-------Fetch nucseq from EBI website-------########
		  
		  use LWP::Simple;
		  @content = get("$ebi_url") or die "\nCan't get content from $ebi_url\n";
		  @split_content = split(/\n/, "@content");
		  $i_max = $#split_content;
		  
		  ##########################################
		  #-------GET TRANSLATION TBL NUMBER-------#
		  $transl_tbl="";
		  foreach(@split_content){
		      if(/\/transl_table=([\d]+)/){
			  $transl_tbl=$1;
		      }
		  }
		  $transl_tbl="1" if($transl_tbl eq ""); #If transl tbl nb is not mentionned, its because transl tbl 1 is used
		  #print "\nTransl_tbl: $transl_tbl";
		  
		  #########################################
		  #------------GET TRANSLATION------------#
		  
		  #Search line at which translation starts
		  for($i=0; $i<=$i_max; $i++){
		      if($split_content[$i] =~ /\/translation="/i) {
			  $translation_line_start = $i;
			  last;
		      }
		  }
		  
		  #First, check if translation starts and stops at the same line
		  $solved = "2";
		  for($i=$translation_line_start; $i<=$i_max; $i++){
		      if($split_content[$i] =~ /\/translation="([\w]+)"/i){
			  $translation = $1;
			  #print ALIGN "\nncbi_aa=$translation"; #Write translation into the file
			  $solved = "1";
			  last;
		      }
		  }
		  
		  #If not, find where translation stops
		  if($solved != "1"){
		      for($i=($translation_line_start + 1); $i<$i_max; $i++){
			  if($split_content[$i] =~ /"/){
			      $translation_line_stop = ($i);
			      last;
			  }
		      }
		      
		      #Store the first line into $translation
		      $position = (index($split_content[$translation_line_start], "\"")) + 1;
		      $firstline = substr($split_content[$translation_line_start], $position);
		      $translation = $firstline;		
		      
		      #Add the middle lines to $translation
		      $i = (($translation_line_start) + 1);
		      while($i < $translation_line_stop){
			  $line = substr($split_content[$i], 2);
			  if($line =~ /\w/i){
			      $translation .= "$&$'";
			  }
			  $i++
		      }
		      
		      #Add the last line to $translation
		      $lastline = substr($split_content[$translation_line_stop], 2);
		      if($lastline =~ /\s([\w]+)"/i){
			  $translation .= "$1";
			  #print "\nTranslation stored.";
		      }	       
		  }
		  
		  #see if the pfam_aa can be matched in the translation if yes get the dna sequence
		  if($translation =~ /($pfam_aa)/){
		      #print "\nFound pfam_aa in translation.";
		      
		      ##########################################
		      #------------GET DNA SEQUENCE------------#
		      
		      #Search line at which DNA sequence starts
		      for($i=0; $i<=$i_max; $i++){
			  if($split_content[$i] =~ /SQ\s+Sequence/i) {
			      $seq_line_start = ($i +1);
			      last;
			  }
		      }
		      
		      #Search line at which DNA sequence stops
		      for($i=$seq_line_start; $i<=$i_max; $i++){
			  if($split_content[$i] =~ /\/\//){
			      $seq_line_stop = ($i);
			      last;
			  }
		      }
		      
		      #Store DNA sequence into variable
		      $i = $seq_line_start;
		      while($i < $seq_line_stop){ #for each line, do the following:
			  $_ = $split_content[$i];
			  s/\s+//g; #remove spaces
			  s/\d+//g; #remove nucleotide position
			  $dna_seq .= $_; #write the line in the file
			  $i++
		      }
		      $dna_seq = uc($dna_seq);
		      
		      #Find the nucleotide sequence in dna_seq that codes for the pfam_aa
		      @coding_seq = ();
		      $nuc_tot = length($dna_seq);
		      $pfam_aa_length = length($pfam_aa);
		      $i = 0; #reset position to 0 in dna_seq
		      $pfam_aa_pos = 0; #reset position to 0 in pfam_aa seq
		      $found_match_pos = 0;
		      while($i < $nuc_tot){ #parse each codon of dna_seq until it finds the seq that matches the PFAM_AA
			  $codon = substr($dna_seq,$i,3); #store codon
			  $aa_from_nuc_seq = $transl_tbl{$codon.$transl_tbl};
			  $aa = substr($pfam_aa,$pfam_aa_pos,1); #store aa
			  if($aa_from_nuc_seq eq $aa){ #check if codon gives same aa as the aa from pfam_aa
			      $found_match_pos = $i if($found_match_pos == 0);
			      push @coding_seq, $codon;
			      $pfam_aa_pos++;
			      if($pfam_aa_pos == $pfam_aa_length){ #stop if we reached the end of the aa seq
				  $nuc_seq_matched=1;
				  last;
			      }
			      $i += 3; #continue on to next codon
			  }
			  else{ #if there's no match: go back to begining of aa seq and clear coding_seq array
			      if($found_match_pos !=0){
				  $i = $found_match_pos+3;
				  $found_match_pos = 0;
			      }
			      else{
				  $i += 3;
			      }
			      $pfam_aa_pos = 0;
			      @coding_seq = ();
			  }
		      }
		      
		      #Store nucleotide sequence in a string
		      foreach(@coding_seq){
			  $coding_seq_string .= $_; 
		      }
		      
		      #Check if we haven't found the nucleotide sequence that codes for the pfam_aa chain
		      if($nuc_seq_matched == 0){ #Error #4
			  #print "\nCouldn't find nucleotide sequence corresponding to aa seq!";
			  $protein_nuc_no_match++;
			  push @protein_nuc_no_match, $uniprot_id;
			  next PROTEIN;
		      }
		      
		      $nb_codons_seq = (length($coding_seq_string))/3;
		      
		      #If nb of codons = nb of aa in the pfam alignment store all info into the file
		      if(($pfam_aa_length == $nb_codons_seq) && $nuc_seq_matched ==1){
			  $protein_stored++;
			  
			  $coding_seq_string = uc($coding_seq_string);
			  #print "\nNb of codons matched nb of aa!\n";
			  
			  #Write data into DATA file
			  print ALIGN "u_id=$uniprot_id//";
			  print ALIGN "\norg_id=$org_id//";
			  print ALIGN "\npfam_ori=$pfam_aa_ori//";
			  print ALIGN "\npfam_aa=$pfam_aa//";
			  print ALIGN "\ntransl_tbl=$transl_tbl//";
			  print ALIGN "\ncoding_seq=$coding_seq_string//\n";
		      }
		      elsif(($pfam_aa_length != $nb_codons_seq) && $nuc_seq_matched ==1){ #Error #5
			  #print "\nNb of codons DIDN'T match nb of aa!";
			  $protein_wrong_codon_nb++;
			  push @protein_wrong_codon_nb, $uniprot_id;
			  next PROTEIN;
		      }
		  }
		  else{ #Error #3
		      #print "\nCounld't find pfam_aa in translation!";
		      $protein_aa_not_matched++;
		      push @protein_aa_not_matched, $uniprot_id;
		      next PROTEIN;
		  }
	      }
	  }
	  $nuc_seq_matched=0;
	  $matched = 0;
	  $translation = "";	  
	  $translation_line_start = "";
	  $translation_line_stop = "";
	  $dna_seq = "";
	  $seq_line_start = "";
	  $seq_line_stop = "";
	  $coding_seq_string = "";
	  $pfam_aa_length = "";
	  $ebi_url = "";
	  $nb_codons_seq = "";
	  $i_max = "";
      }
	close CONTENT;
	close ALIGN;
	
	###------DATA EXTRACTION STATISTICS------###
	#print "\n\n------------------------------------";
	#print "\n-----DATA EXTRACTION STATISTICS-----";
	#print "\n------------------------------------";
	print "$protein_stored proteins over $protein_count passed stage 1. See extraction log file for details.";
	
	if($protein_stored > 0){
	    open LOG_ERRORS, ">results/$store_res_dir/extraction_log" or die "Cant open error log file";
	    print LOG_ERRORS "Proteins with no codon usage info\n" if($#protein_no_cu>=0);
	    print LOG_ERRORS "$_ $pfam_file\n" foreach(@protein_no_cu);
	    
	    print LOG_ERRORS "\nProteins for which the nucleotide sequence couldn't be retrived\n" if($#protein_no_ebi>=0);
	    print LOG_ERRORS "$_ $pfam_file\n" foreach(@protein_no_ebi);
	    
	    print LOG_ERRORS "\nProteins for which the aa sequence in the pfam alignment couldn't be matched with the one on the uniprot website\n" if($#protein_aa_not_matched>=0);
	    print LOG_ERRORS "$_ $pfam_file\n" foreach(@protein_aa_not_matched);
	    
	    print LOG_ERRORS "\nProteins with wrong number of codons\n" if($#protein_wrong_codon_nb>=0);
	    print LOG_ERRORS "$_ $pfam_file\n" foreach(@protein_wrong_codon_nb);
	    
	    print LOG_ERRORS "\nProteins for which the aa seq from the alignement wasn't matched with the nucleotide sequence\n" if($#protein_nuc_no_match>=0);
	    print LOG_ERRORS "$_ $pfam_file\n" foreach(@protein_nuc_no_match);
	    
	    close LOG_ERRORS;
	}
	
	$file_path = "stage1/pfam_".$pfam_file;
	unlink($file_path) if($protein_stored == 0);
	
	@pfam_file=();
	@protein_no_cu=();
	@protein_no_ebi=();
	@protein_aa_not_matched=();
	@protein_wrong_codon_nb=();
	@protein_nuc_no_match=();
	$protein_stored=0;
	
    }
    
    print "\n\n--------------- STAGE 2 ---------------\n";
    
    $max_cut_index=$#cut;
    for($i=0; $i<=$max_cut_index; $i++){
	if($cut[$i] =~ /org_id=([0-9]+)\|/){
	    $org_id=$1;
	    $j=$i+1;
	    $cu_hash{$org_id}=$cut[$j];
	}
    }
    
    @pfam_output_list_glob = glob "stage1/*";
    print "Assigning codon usage frequencies to every codon\n";
  FILE: foreach(@pfam_output_list_glob){ #foreach alignment
      $pfam_file_path = $_;
      $pfam_file_name = $1 if($pfam_file_path =~ /\/pfam_([A-Z_0-9\.]+)$/i);
      print "\n$pfam_file_name\n";
      open OUTPUT, ">>stage2/$pfam_file_name" or die "cant open file to freq";
      open DATA, "<$pfam_file_path" or die "\nCan't open data file\n";
      #########################################
      ##########   PARSE PFAM FILE   ##########
      #########################################
    LINE: foreach(<DATA>){
	if($_ =~ /u_id=(.+)\/\//){
	    $u_id = $1; #Store u_id
	    next LINE;
	}
	
	if($_ =~ /org_id=([0-9]+)\/\//){
	    $org_id = $1; #Store the organism's id
	    next LINE;
	}
	
	if($_ =~ /pfam_ori=(.*)\/\//){ #Store line from the alignment that corresponds to this protein
	    $pfam_ori = $1;
	    $pfam_ori_length = length($pfam_ori);
	    next LINE;
	}
	
	if($_ =~ /transl_tbl=(.*)\/\//){ #Store the translation table number
	    $transl_tbl = $1;
	    next LINE;
	}
	
	if($_ =~ /coding_seq=(.*)\/\//){ #Get nucleotide sequence and start rc analysis
	    $coding_seq = $1;
	    
	    %codon_freq=('TTT'=>'NA', 'TTC'=>'NA', 'TTA'=>'NA', 'TTG'=>'NA', 'TCT'=>'NA', 'TCC'=>'NA', 'TCA'=>'NA', 'TCG'=>'NA', 'TAT'=>'NA', 'TAC'=>'NA', 'TAA'=>'NA', 'TAG'=>'NA', 'TGT'=>'NA', 'TGC'=>'NA', 'TGA'=>'NA', 'TGG'=>'NA', 'CTT'=>'NA', 'CTC'=>'NA', 'CTA'=>'NA', 'CTG'=>'NA', 'CCT'=>'NA', 'CCC'=>'NA', 'CCA'=>'NA', 'CCG'=>'NA', 'CAT'=>'NA', 'CAC'=>'NA', 'CAA'=>'NA', 'CAG'=>'NA', 'CGT'=>'NA', 'CGC'=>'NA', 'CGA'=>'NA', 'CGG'=>'NA', 'ATT'=>'NA', 'ATC'=>'NA', 'ATA'=>'NA', 'ATG'=>'NA', 'ACT'=>'NA', 'ACC'=>'NA', 'ACA'=>'NA', 'ACG'=>'NA', 'AAT'=>'NA', 'AAC'=>'NA', 'AAA'=>'NA', 'AAG'=>'NA', 'AGT'=>'NA', 'AGC'=>'NA', 'AGA'=>'NA', 'AGG'=>'NA', 'GTT'=>'NA', 'GTC'=>'NA', 'GTA'=>'NA', 'GTG'=>'NA', 'GCT'=>'NA', 'GCC'=>'NA', 'GCA'=>'NA', 'GCG'=>'NA', 'GAT'=>'NA', 'GAC'=>'NA', 'GAA'=>'NA', 'GAG'=>'NA', 'GGT'=>'NA', 'GGC'=>'NA', 'GGA'=>'NA', 'GGG'=>'NA');
	    
	    @split_cu=();
	    @split_cu = split(/\|/, $cu_hash{$org_id});
	    foreach(@split_cu){
		$codon_freq{$1} = $2 if($_ =~ /^([a-z]{3})([0-9\.]+)/i);
	    }
	    
	    print OUTPUT "$u_id,$org_id:";
	    ###################################################
	    ###---------------------------------------------###
	    #### For every column position of this protein ####
	    ###---------------------------------------------###
	    ###################################################
	    $count = 0; #keeps track of the aa position in the alignment
	    for($i=0; $i<$pfam_ori_length; $i++){ #Here $i represents the column number of the alignment
		$pfam_char = substr($pfam_ori,$i,1);
		if($pfam_char =~ /[a-z]{1}/i){ #If we encounter an aa in the pfam alignment
		    $count++;
		    $nuc_seq_pos = (($count - 1)*(3)); #position of corresponding codon in coding_seq
		    $codon = substr($coding_seq,$nuc_seq_pos,3); #Fetch the codon
		    $aa=$transl_tbl{$codon.$transl_tbl};
		    if($codon_freq{$codon} ne "NA"){
			print OUTPUT "$aa,$codon_freq{$codon},$codon-";
		    }
		    else{
			print OUTPUT "x-";
		    }
		}
		else{
		    print OUTPUT "x-";
		}
	    }
	    $org_id="";
	    $u_id="";
	    $pfam_ori="";
	    $transl_tbl="";
	    $coding_seq="";
	    print OUTPUT "\n";
	}
	
    }#done for this protein
      close OUTPUT;
  }#done for this pfam
}

print "\n\n--------------- STAGE 3 ---------------\n";
print "Analyse for rare codon clusters and predict unstructured regions\n";

$pfam_w_pause=0;
$pfam_w_cluster=0;
$tot_nb_pause=0;
$tot_nb_cluster=0;

@avrg_files_glob = glob "stage2/*";

#Get nb of positions at the left and right of the window (left is smaller than right if unequal window)
$half_win=($win_size-1)/2;
if($half_win=~/\./){
    $right=$1 if($half_win=~/([0-9]+)\./);
    $left=($win_size-1)-$right;
}else{
    $right=$left=$half_win;
}

mkdir "results/$store_res_dir/html" or die "cant create html dir";
system("cp files/style.css results/$store_res_dir/html/style.css");
system("cp files/sherlocc.jpg results/$store_res_dir/html/sherlocc.jpg");

####-----------------------------------------------####
######-----    START PARSING PFAM FILES     -----######
####-----------------------------------------------####

#Parse each file in the given directory
PFAM: foreach(@avrg_files_glob){
    @col_freq=();
    @pauses=();
    $pfam_pauses="";
    @col_type=();
    @fusionned_clusters=();

    $cluster_positions="";
    $this_pfam_pause=0; #will change to 1 if a pause is encountered
    $this_pfam_cluster=0; #will change if a cluster is found
    $pfam_file_path = $_;
    if($pfam_file_path =~ /\/([a-z0-9\.]+)$/i){
	$pfam_file_name = $1;
	$pfam_name = $1;
    }
    $pfam_name=~s/\.(.+)$//;
    $html_file=$pfam_name.".html";
    
    #open file to store results
    open RESULTS, ">>results/$store_res_dir/results" or die "cant create sherlocc result file";

    $alignment_length=0;
    #Get alignment length
    open DATA, "<$pfam_file_path" or die "\nCan't open data file\n";
    foreach(<DATA>){
	$_ =~ s/-$//;
	if($_ =~ /^(.+):([-\.a-z,0-9]+)$/i){
	    @array=split(/-/, $2);
	    $alignment_length=$#array+1;
	}
	last;
    }
    close DATA;

    $prot_nb_in_align=0;
    #Get nb of proteins
    open DATA, "<$pfam_file_path" or die "\nCan't open data file\n";
    foreach(<DATA>){
	$prot_nb_in_align++;
    }
    close DATA;

    #Create and start printing HMTL file
    if($cr_html eq "y"){
	open HTML, ">>results/$store_res_dir/html/$html_file" or die "cant open html output";
	print HTML "<html><head><title>$pfam_file_name</title><link rel=stylesheet type=text/css href=style.css /></head>";
	print HTML "<body topmargin=10 bottommargin=10 leftmargin=10 rightmargin=10><img src=sherlocc.jpg><br><br>";
	print HTML "<span class=link><a target=_new href=http://pfam.sanger.ac.uk/family/$pfam_name>$pfam_file_name</a></span>";
	print HTML "<br><br><span class=corps>codon usage average threshold: $avrg_lim";
	print HTML "<br>window size for unstructured/structured predictions: $win_size";
	print HTML "<br>Number of proteins in this alignment: $prot_nb_in_align";
	print HTML "<br>Residue length of the alignment: $alignment_length</span><br>";
	print HTML "<br><br><table border=0 cellpadding=4 cellspacing=0>";
	
	#Print first row of HTML alignement table (row with column positions)
	open DATA, "<$pfam_file_path" or die "\nCan't open data file\n";
	foreach(<DATA>){
	    $_ =~ s/-$//;
	    if($_ =~ /^(.+):([-\.a-z,0-9]+)$/i){
		@array=split(/-/, $2);
		print HTML "<tr>";
		print HTML "<td align=center>&nbsp;</td>";
		for($i=0;$i<=$#array;$i++){
		    $pos=$i+1;
		    print HTML "<td align=center><font size=1><u><b>$pos</b></u></font></td>";
		}
		print HTML "</tr>";
		last;
	    }
	}
	close DATA;
    }
    
    #Build the matrix and print the other HTML rows
    $row=0;
    open DATA, "<$pfam_file_path" or die "\nCan't open data file\n";
    foreach(<DATA>){
	@protein_seq=();
	
	$_ =~ s/-$//;
	if($_ =~ /^(.+),(.+):([-\.a-z,0-9]+)$/i){
	    $u_id=$1;
	    $org_id=$2;
	    @array=split(/-/, $3);
	    
	    #create 1D protein sequence matrix on which to run findex algorithm
	    $col=0;
	    foreach(@array){ #foreach position of this alignement row
		if($_ =~ /([a-z]{1}),([\.0-9]+),([a-z]+)/i){
		    $this_aa=uc($1);
		    $protein_seq[$col]=$this_aa;
		}elsif($_ =~ /x/i){
		    $protein_seq[$col]="x";
		}
		$col++;
	    }
	    
	    #print html row with aa and gaps
	    print HTML "<tr>" if($cr_html eq "y");
	    print HTML "<td align=center bgcolor=E8E8E7><font size=1><a target=_new href=http://www.uniprot.org/uniprot/$u_id>$u_id</a><br><a target=_new href=http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=$org_id>$org_id</a></font></td>" if($cr_html eq "y");
	    $col=0;
	    foreach(@array){
		if($_ =~ /([a-z]{1}),([\.0-9]+),([a-z]+)/i){
		    $this_aa=uc($1);
		    print HTML "<td bgcolor=#E1E1E1 align=center><font size=1>$1 <b>$2</b>&nbsp;$3</font></td>" if($cr_html eq "y");
		    $matrix[$row][$col]=$2;
		}
		elsif($_ =~ /x/){
		    print HTML "<td bgcolor=#E1E1E1 align=center><font size=1>x</font></td>" if($cr_html eq "y");
		    $matrix[$row][$col]="x";
		}
		$col++;
	    }
	    print HTML "</tr>" if($cr_html eq "y");
	    
	    #||||||-------------------------------------------------------------------|||||||
	    #||||||-------------------------Foldindex algo----------------------------|||||||
	    #||||||-------------------------------------------------------------------|||||||

	    for($x=0; $x<$alignment_length; $x++){
		if(($x>=$left) && ($x<($alignment_length-$right))){
		    $hydro_sum=0;
		    $charge_sum=0;
		    $mean_hydro=0;
		    $mean_charge=0;
		    $nb_gaps_in_this_win=0;
		    $value="";
		    
		    #Sum left part of the window
		    for($l=$x-$left; $l<$x; $l++){
			if($protein_seq[$l] eq "x"){
			    $nb_gaps_in_this_win++;
			}else{
			    $hydro_sum+=$hydro{$protein_seq[$l]};
			    $charge_sum+=$charge{$protein_seq[$l]};
			}
		    }
		    
		    #Sum center position of the window
		    if($protein_seq[$x] eq "x"){
			$nb_gaps_in_this_win++;
		    }else{
			$hydro_sum+=$hydro{$protein_seq[$x]};
			$charge_sum+=$charge{$protein_seq[$x]};
		    }
		    
		    #Sum right part of the window
		    for($r=$x+1; $r<=($x+$right); $r++){
			if($protein_seq[$r] eq "x"){
			    $nb_gaps_in_this_win++;
			}else{
			    $hydro_sum+=$hydro{$protein_seq[$r]};
			    $charge_sum+=$charge{$protein_seq[$r]};
			}
		    }
		    
		    $mean_hydro=$hydro_sum/($win_size-$nb_gaps_in_this_win) unless(($win_size-$nb_gaps_in_this_win) <= 0);
		    $charge_sum=~s/-//;
		    $mean_charge=$charge_sum/($win_size-$nb_gaps_in_this_win) unless(($win_size-$nb_gaps_in_this_win) <= 0);
		    $mean_hydro=sprintf("%.3f", $mean_hydro);
		    $mean_charge=sprintf("%.3f", $mean_charge);
		    
		    #Get window interval
		    $left_int=$x-$left+1;
		    $right_int=$x+$right+1;
		    
		    if(($win_size-$nb_gaps_in_this_win) > 0){
			$value=(2.785*($mean_hydro))-($mean_charge)-(1.151);
			$value=sprintf("%.3f", $value);
			$findex_matrix[$row][$x]=$value." ".$left_int."-".$right_int.", ".$nb_gaps_in_this_win;
		    }else{
			$findex_matrix[$row][$x]="na";
		    }
		}else{
		    $findex_matrix[$row][$x]="na";
		}
	    }
	}
	print HTML "<tr>" if($cr_html eq "y");
	print HTML "<td align=center bgcolor=E8E8E7><font size=1></font></td>" if($cr_html eq "y");
	#------------------------------------
	#Print fdx row for each protein
	for($x=0; $x<($alignment_length); $x++){
	    $bgcolor="#FFFFFF";
	    if($findex_matrix[$row][$x] eq "na"){
		print HTML "<td align=center width=35 bgcolor=#FFFFFF><font size=1>na</font></td>" if($cr_html eq "y");
	    }else{
		if($findex_matrix[$row][$x]=~/^([\.0-9-]+\s+)/){
		    $findex=$1;
		    if($findex=~/-/){
			$bgcolor="efa8a8";
		    }else{
			$bgcolor="#d0f2bb";
		    }
		}
		print HTML "<td align=center width=35 bgcolor=$bgcolor><font size=1>$findex</font></td>" if($cr_html eq "y");
	    }
	}
	print HTML "</tr>" if($cr_html eq "y");

	$row++;
    }
    close DATA; #close get_avrg file

    $tresh=($row*(7))*(0.75);

    #########################    SHERLOCC    #############################
    ######################################################################
    ###      ...|        CALCULATE 7 CODON-WIDE WINDOW        |...     ###
    ######################################################################
    print HTML "<tr>" if($cr_html eq "y");
    print HTML "<td align=center bgcolor=C7C7C6><font size=1><b>Average</b></font></td>" if($cr_html eq "y");
    for($i=0; $i<$col; $i++){ #loop every position
	$freq_sum=0;
	$freq_count=0;
	$na_count=0;
	for($r=0; $r<$row; $r++){ #calculate frequency avrg for a 7 codon-wide window
	    COL: for($c=($i-3); $c<=($i+3); $c++){
		if(($c < 0) || ($c > (($col)-1))){
		    $na_count++;
		    next COL;
		}
		if($matrix[$r][$c] =~ /[0-9\.]+/){
		    $freq_sum+=$matrix[$r][$c];
		    $freq_count++;
		}
		elsif($matrix[$r][$c] =~ /x/){
		    $na_count++;
		}
	    }
	}
	if($freq_count >= $tresh){ #we have a representative col
	    $avrg=$freq_sum/$freq_count;
	    $col_freq[$i]=$avrg;
	    $rounded_avrg=sprintf("%.2f", $avrg);
	    if($avrg <= $avrg_lim){ #We have a pause
		$col_type[$i]="p";
		$tot_nb_pause++;
		$this_pfam_pause=1;
		print HTML "<td bgcolor=#C45700 align=center><font size=1><b>$rounded_avrg</b></font></td>" if($cr_html eq "y");
		push @pauses, $i;
	    }
	    else{ #we dont have a pause
		$col_type[$i]="c";
		print HTML "<td bgcolor=C7C7C6 align=center><font size=1><b>$rounded_avrg</b></font></td>" if($cr_html eq "y");
	    }
	}
	else{ #col is not representative
	    $col_type[$i]="c";
	    print HTML "<td bgcolor=C7C7C6 align=center><font size=1>NA</font></td>" if($cr_html eq "y");
	}
    }
    print HTML "</tr>" if($cr_html eq "y");

    foreach(@pauses){
	$_+=1;
	$pfam_pauses.="$_,";
    }
    $pfam_pauses=~s/,$//;

    ############################    SHERLOCC    ##############################
    ##         .|#################################################|.        ##
    ##       ...|########---LOCATE AND COMBINE CLUSTERS---########|...      ##
    ##     .....|#################################################|.....    ##
    ##########################################################################

    if($this_pfam_pause==1){ #we have at least one pause -> lets look for clusters...
	@clusters=();
	$pfam_w_pause++;
	#Parse for significant windows (clusters)
	for($i=0; $i<=($#col_type-6); $i++){
	    $win_pause_count=0;
	    $win_pauses="";

	    for($j=$i; $j<=($i+6); $j++){ #count nb of pauses for this window
		if($col_type[$j] eq "p"){
		    $win_pause_count++;
		    $win_pauses.="$j,";
		}
	    }
	    if($win_pause_count >= 4){ #we found a significant cluster
		$this_pfam_cluster=1;
		@split_pauses=split(/,/, $win_pauses);
		$start=$split_pauses[0];
		$stop=$split_pauses[$#split_pauses];
		push @clusters, "$start-$stop";
	    }
	}
	
	if($this_pfam_cluster==1){ #We have at least one cluster
	    $this_pfam_cluster_freq_sum=0; #set cluster_frequency_sum of this pfam to zero
	    $this_pfam_cluster_freq_avrg=0; #set cluster_frequency_avrg of this pfam to zero
	    $this_pfam_cluster_count=0;
	    print RESULTS "pfam=$pfam_file_name//\nclusters=";
	    $pfam_w_cluster++;
	    $nb_ind_clus=0;
	    foreach(@clusters){$nb_ind_clus++;} #count nb of clusters

	    ##############----  WE HAVE ONLY 1 CLUSTER  ----############
	    if($nb_ind_clus==1){
		$tot_nb_cluster++;
		if($clusters[0] =~ /([0-9]+)-([0-9]+)/){
		    $start=$1;
		    $stop=$2;
		    $pos_count=0;
		    $freq_sum=0;
		    $cluster_avrg=0;
		    for($k=$start; $k<=$stop; $k++){
			$freq_sum+=$col_freq[$k];
			$pos_count++;
		    }
		    $cluster_avrg=$freq_sum/$pos_count;
		    $this_pfam_cluster_freq_sum+=$cluster_avrg;
		    $this_pfam_cluster_count++;
		}
		$start_adj=($start+1);
		$stop_adj=($stop+1);
		print RESULTS "$start_adj-$stop_adj|$cluster_avrg";
		push @fusionned_clusters, "$start_adj-$stop_adj|$cluster_avrg";
		for($cluster_col=($start_adj-1); $cluster_col<=($stop_adj-1); $cluster_col++){
		    $cluster_positions.=",$cluster_col,";
		}
		$this_pfam_cluster_freq_avrg+=$this_pfam_cluster_freq_sum/$this_pfam_cluster_count;
		print RESULTS "//\nnb_of_rare_codon_cluster=$this_pfam_cluster_count//\nclusters_freq_avrg=$this_pfam_cluster_freq_avrg//";
	    }
	    ##############----  WE HAVE MORE THAN 1 CLUSTER  ----############
	    elsif($nb_ind_clus>1){
		if($clusters[0] =~ /([0-9]+)-([0-9]+)/){ #get positions of first cluster
		    $start=$1;
		    $stop=$2;
		}
		$last_cluster_index=$#clusters; #get array index of last cluster
		for($i=1; $i<=$last_cluster_index; $i++){
		    if($clusters[$i] =~ /([0-9]+)-([0-9]+)/){ #get positions of this cluster
			$this_start=$1;
			$this_stop=$2;
		    }
		    if(($this_start >= $start) && ($this_start <= $stop)){ #if they overlap define new stop position
			$stop=$this_stop;
		    }
		    elsif($this_start>$stop){ #if they dont overlap store the first one and define new start/stop positions
			$cluster_int=$start."-".$stop;
			if($cluster_int =~ /([0-9]+)-([0-9]+)/){
			    $start=$1;
			    $stop=$2;
			    $pos_count=0;
			    $freq_sum=0;
			    $cluster_avrg=0;
			    for($k=$start; $k<=$stop; $k++){
				$freq_sum+=$col_freq[$k];
				$pos_count++;
			    }
			    $cluster_avrg=$freq_sum/$pos_count;
			    $this_pfam_cluster_freq_sum+=$cluster_avrg;
			    $this_pfam_cluster_count++;
			}
			$start_adj=($start+1);
			$stop_adj=($stop+1);
			print RESULTS "$start_adj-$stop_adj|$cluster_avrg,";
			push @fusionned_clusters, "$start_adj-$stop_adj|$cluster_avrg";
			for($cluster_col=($start_adj-1); $cluster_col<=($stop_adj-1); $cluster_col++){
			    $cluster_positions.=",$cluster_col,";
			}
			$tot_nb_cluster++;
			#new reference positions
			$start=$this_start;
			$stop=$this_stop;
		    }
		    
		    if($i==$last_cluster_index){ #store the cluster if its the last one
			$cluster_int=$start."-".$stop;
			if($cluster_int =~ /([0-9]+)-([0-9]+)/){
			    $start=$1;
			    $stop=$2;
			    $pos_count=0;
			    $freq_sum=0;
			    $cluster_avrg=0;
			    for($k=$start; $k<=$stop; $k++){
				$freq_sum+=$col_freq[$k];
				$pos_count++;
			    }
			    $cluster_avrg=$freq_sum/$pos_count;
			    $this_pfam_cluster_freq_sum+=$cluster_avrg;
			    $this_pfam_cluster_count++;
			}
			$start_adj=($start+1);
			$stop_adj=($stop+1);
			print RESULTS "$start_adj-$stop_adj|$cluster_avrg,";
			push @fusionned_clusters, "$start_adj-$stop_adj|$cluster_avrg";
			for($cluster_col=($start_adj-1); $cluster_col<=($stop_adj-1); $cluster_col++){
			    $cluster_positions.=",$cluster_col,";
			}
			$tot_nb_cluster++;
		    }
		}
		$this_pfam_cluster_freq_avrg+=$this_pfam_cluster_freq_sum/$this_pfam_cluster_count;
		print RESULTS "//\nnb_of_rare_codon_cluster=$this_pfam_cluster_count//\nclusters_freq_avrg=$this_pfam_cluster_freq_avrg//\n";
	    }
	}	
    } #Finished combining and fusionning clusters

    #------------------------------------------------------------------------
    #PRINT CLUSTER BAR
    if($cr_html eq "y"){
	print HTML "<tr>";
	print HTML "<td height=10 align=center bgcolor=ffffff></td>";
	for($i=0; $i<$col; $i++){
	    print HTML "<td height=10 bgcolor=ffffff></td>";
	}
	print HTML "</tr>";
	print HTML "<tr>";
	print HTML "<td align=center bgcolor=949493 nowrap><font size=1><b>Rare Codon Clusters</b></font></td>";
	for($i=0; $i<$col; $i++){
	    if($cluster_positions =~ /,$i,/){
		print HTML "<td bgcolor=#A80302></td>";
	    }else{
		print HTML "<td bgcolor=3cb821>&nbsp;</td>";
	    }
	}
	print HTML "</tr>";
    }

    #Print HTML spacer before findex avrg row
    if($cr_html eq "y"){
	print HTML "<tr><td height=10 align=center bgcolor=#FFFFFF></td>";
	for($x=0; $x<$alignment_length; $x++){print HTML "<td bgcolor=#FFFFFF></td>";}
	print HTML "</tr>";
    }

    #Print HTML row of RAW FINDEX
    print HTML "<tr>" if($cr_html eq "y");
    print HTML "<td align=center bgcolor=949493 nowrap><font size=1><b>Unstructured index values</b></font></td>" if($cr_html eq "y");
    for($x=0; $x<$alignment_length; $x++){
	#Get fdx_value for this position
	if(($x>=$left) && ($x<($alignment_length-$right))){ #If we have a value for this column (inside window range)
	    $col_findex_avrg=0;
	    $col_sum=0;
	    $col_count=0;
	    for($y=0; $y<$row; $y++){ #Average the value over the column
		if($findex_matrix[$y][$x]=~/^([-0-9\.]+)\s+/){
		    $col_sum+=$1;
		    $col_count++;
		}
	    }
	    $col_findex_avrg=$col_sum/$col_count if($col_count>0);
	    $col_findex_avrg=sprintf("%.3f", $col_findex_avrg);
	}else{
	    $col_findex_avrg="";
	}
	
	#Determine cell bg_color
	if($col_findex_avrg=~/-/){
	    $bgcolor="efa8a8"; #below 0 (unfolded)
	}elsif($col_findex_avrg eq ""){
	    $bgcolor="#DFDFDF"; #no value
	}else{
	    $bgcolor="#d0f2bb"; #positive (folded)
	}

	print HTML "<td bgcolor=$bgcolor><font size=1>$col_findex_avrg</font></td>" if($cr_html eq "y");
	$fx_col[$x]=$col_findex_avrg;
    }
    print HTML "</tr>" if($cr_html eq "y");

    #Print HTML spacer before findex avrg row
    if($cr_html eq "y"){
	print HTML "<tr><td height=10 align=center bgcolor=#FFFFFF></td>";
	for($x=0; $x<$alignment_length; $x++){print HTML "<td bgcolor=#FFFFFF></td>";}
	print HTML "</tr></table>";
    }

    #########################################     FINDEX      ###########################################
    ##    ...|########---LOCATE AND COMBINE CLUSTERS of raw and cons UNFOLDED REGIONS---########|...   ##
    #####################################################################################################

    @fx_ints=();
    @cons_fx_ints=();
    @mp_unf=();
    @fx_conservative=();
    $nb_cons_ints=0;
    $nb_fx_ints=0;
    for($x=0; $x<=$alignment_length; $x++){
	$fx_conservative[$x]="f";
    }
    $count_unf=0;
    FX: for($x=$left; ($x<($alignment_length-$right)); $x++){
	next FX if($fx_col[$x] eq "");
	if(($fx_col[$x]<0) && ($count_unf==0)){ #if we encounter a first cons position
	    $count_unf++;
	    $start_pos=$x;
	}elsif(($fx_col[$x]<0) && ($count_unf>0) && ($x<($alignment_length-$right-1))){ #if we encounter other cons pos
	    $count_unf++;
	}elsif(($count_unf>0) && ($x==($alignment_length-$right-1))){ #if its the end of the alignment and we have a clus opened
	    $end=$x+1;
	    $start_pos+=1;
	    push @fx_ints, $start_pos."-".$end;
	    $nb_fx_ints++;
	    $width=$end-$start_pos+1;
	    
	    if($width>25){
		$nb_cons_ints++;
		push @cons_fx_ints, $start_pos."-".$end; #push in conservative interval array

		for($x=$start_pos-1; $x<=$end-1; $x++){
		    $fx_conservative[$x]="unf"; #update conservative unfolded position array
		}

		$middle_point=$end-(($end-$start_pos+1)/2); #Get middle point
		push @mp_unf, $middle_point; #Store middle point of this cons unfolded interval
	    }
	    
	    undef($end);
	    $count_unf=0;
	}elsif(($fx_col[$x]>0) && ($count_unf>0)){
	    $end=$x;
	    $start_pos+=1;
	    push @fx_ints, $start_pos."-".$end; #push in raw_interval array
	    $width=$end-$start_pos+1;
	    $nb_fx_ints++;
	    
	    if($width>25){ #its a conservative interval
		$nb_cons_ints++;
		push @cons_fx_ints, $start_pos."-".$end; #push in conservative interval array

		for($x=$start_pos-1; $x<=$end-1; $x++){
		    $fx_conservative[$x]="unf"; #update conservative unfolded position array
		}

		$middle_point=$end-(($end-$start_pos+1)/2); #Get middle point
		push @mp_unf, $middle_point; #Store middle point of this cons unfolded interval
	    }
	    
	    undef($end);
	    $count_unf=0;
	}
    }
    undef($middle_point);

    $nb_of_cons_fx_ints=0;
    foreach(@cons_fx_ints){
	$nb_of_cons_fx_ints++;
    }

    $count_clusters=0;
    foreach(@fusionned_clusters){
	$count_clusters++;
    }
    @mp_clus=();
    $total_cluster_length=0;
    #------------------------------------------------------------------------
    #PRINT TABLE WITH CLUSTER POSITIONS AND FREQ AVERAGE
    if($cr_html eq "y"){
	print HTML "<br><br><table border=0 cellpadding=4 cellspacing=0>";
	print HTML "<tr><td width=300 colspan=3 bgcolor=949493><span class=titre>Rare Codon Clusters</span></td></tr>";
	print HTML "<tr><td width=80 bgcolor=#C7C7C6 align=left><span class=corps><b>Position</b></span></td>";
	print HTML "<td width=220 bgcolor=#C7C7C6 align=right><span class=corps><b>Usage Frequency Average</b></span></td>";
	print HTML "<td width=100 bgcolor=#C7C7C6 align=right><span class=corps><b>Middle Point</b></span></td></tr>";
    }
    if(($count_clusters==0) && ($cr_html eq "y")){
	print HTML "<tr><td width=80 align=left><span class=corps>No clusters.</span></td>";
	print HTML "<td width=220 align=right></td>";
	print HTML "<td width=100 align=right></td></tr>";
    }elsif($count_clusters>0){
	foreach(@fusionned_clusters){
	    if($_ =~ /([0-9]+)-([0-9]+)\|([0-9\.]+)/){
		$clus_start=$1;
		$clus_stop=$2;
		$total_cluster_length+=($2-$1+1);
		$avrg = sprintf("%.3f", $3);
		$middle_point=$2-(($2-$1+1)/2);
		if($middle_point=~/\./){
		    $middle_point=~s/\.[0-9]+//g;
		}
		push @mp_clus, $middle_point;
		if($cr_html eq "y"){
		    print HTML "<tr><td width=80 align=left><span class=corps>$clus_start - $clus_stop</span></td>";
		    print HTML "<td width=220 align=right><span class=corps>$avrg</span></td>";
		    print HTML "<td width=100 align=right><span class=corps>$middle_point</span></td></tr>";
		}
		$count_clusters++;
		$last_mp_clus=$middle_point;
	    }
	}
    }

    #Calculate Total CLUS FRACTION
    $clus_fraction=$total_cluster_length/$alignment_length if($alignment_length>0);
    $clus_fraction=sprintf("%.10f", $clus_fraction);
    
    if($cr_html eq "y"){
	print HTML "</table>";
	if($count_clusters>0){
	    print HTML "<p><span class=corps>Total cluster length: $total_cluster_length";
	    print HTML "<br>Middle point of last cluster: $last_mp_clus";
	    print HTML "<br>Fraction of the pfam occupied by rare codon clusters: $clus_fraction</span><br><br>";
	}
    }

    #------------------------------------------------------------------------
    #print unstructured regions
    if($cr_html eq "y"){
	print HTML "<br><br><table border=0 cellpadding=4 cellspacing=0>";
	print HTML "<tr><td width=300 colspan=2 bgcolor=949493><span class=titre>Unstructured regions</span></td></tr>";
	print HTML "<tr><td width=200 bgcolor=#C7C7C6 align=left><span class=corps><b>Position</b></span></td>";
	print HTML "<td width=100 bgcolor=#C7C7C6 align=right><span class=corps><b>Length</b></span></td></tr>";
	if($nb_fx_ints>0){
	    foreach(@fx_ints){
		if($_=~/([0-9]+)-([0-9]+)/){
		    $width=$2-$1+1;
		    if($2!=$1){
			print HTML "<tr><td width=200 align=left><span class=corps>$_</span></td>";
			print HTML "<td width=100 align=right><span class=corps>$width</span></td></tr>";
		    }else{
			print HTML "<tr><td width=200 align=left><span class=corps>$2</span></td>";
			print HTML "<td width=100 align=right><span class=corps>$width</span></td></tr>";
		    }
		}
	    }
	}else{
	    print HTML "<tr><td colspan=2><span class=corps>There is no unstructured region.</span></td></tr>";	
	}
	print HTML "</table>";
    }

    close RESULTS;
    print HTML "</body></html>" if($cr_html eq "y");
    close HTML if($cr_html eq "y");

} # --> Go to next pfam

#Write final statistics into RESULTS FILE
open RESULTS, ">>results/$store_res_dir/results" or die "cant open results file";
print RESULTS "\npfam_with_pauses=$pfam_w_pause\ntotal_nb_of_pauses=$tot_nb_pause\npfam_with_clusters=$pfam_w_cluster\ntot_nb_clusters=$tot_nb_cluster";
close RESULTS;

if($cleanup_dirs eq "y"){
    cleanup("stage1") if(-d "stage1");
    cleanup("stage2") if(-d "stage2");
}

print "\nDone. See results/$store_res_dir/\n\n";

sub cleanup {
    my $dir = shift;
    local *DIR;
    
    opendir DIR, $dir or die "opendir $dir: $!";
    for (readdir DIR) {
	next if /^\.{1,2}$/;
	my $path = "$dir/$_";
	unlink $path if -f $path;
	cleanup($path) if -d $path;
    }
    closedir DIR;
    rmdir $dir or print "error - $!";
}
