#!/usr/bin/perl
 use strict;
 use warnings;


my $ATG_distance = 15;    # number of codons excluded from start codon
my $stop_distance = 10;   # number of codons excluded from stop codon
my $site = "A";   # choose E, P, A and A+1 to calculate the occurrance at different position

#get the CAI reference

my %cai;
 open (hand1, "CAI_AA.txt") or die $!;
 while (<hand1>)   {
 $_=~ s/\s+$// ;
 my ($aa,$co,$rscu,$ca)=split /\t/;
 $cai{$co}="$aa\t$rscu\t$ca";
                   }
 close hand1;

mkdir "speed";

&speed ('sample_name1');
&speed ('sample_name2');
&speed ('sample_name3');
&speed ('sample_name4');

sub speed                                                          {

my $sample = shift;

print "now processing $sample\n\n";

 my $ge; my %gene;    #the ncu number as key and ORF sequence as input
 my %rpf;                    #coverage hash for RPF with ncu and start position as key and score as value
 my %si;                      #size of each gene
 my %rpf_total;           # total number of inframe reads for each gene

    open (hand1, "100_reference_100.fa") or die $!;
    while (<hand1>) {
    $_=~ s/\s+$//;

    if (/^>/)       {
    $ge=$_;
    $ge=~ s/^>//;
    next;}

    $gene{$ge}=$_; }

    close hand1;

 #############
#get the genome length
##############
    foreach my $j (keys %gene)  {
    $si{$j}=length($gene{$j});
                                }

my $start = 0;

if ($site eq "A")  {
	$start = 0; }
elsif ($site eq "P")  {
	$start = -3;  }
elsif ($site eq "E")  {
	$start = -6;}
elsif ($site eq "A+1")  {
	$start = 3;}

######
#get data from RPF files and bed file
#################

  open (hand2, "$sample/$sample\_A.bg") or die $!;

my %rpf_cod_num;  # value of effective codons

  while (<hand2>)                    {
  chomp;
  my ($chr, $st, $end, $value) =split /\t/;
  next if $value == 0;
  for my $i ($st..$end-1)  {
    my $pos = $i + $start;
    $rpf{$chr."\t".$pos}=$value;
    # recording the total reads that are within the coding region and in frame
    if ( $st > 99 + $ATG_distance*3 + $start and $end < ($si{$chr} - 100 - $stop_distance*3 - $start) ) {
	  if (  ($st-100)%3 == 0  )   {
		$rpf_total{$chr} += $value;
		$rpf_cod_num{$chr} ++; }
      }     }40
    }
  close hand2;

#calculate the average RPF density of each gene
my %rpf_avg_gene_dens;

foreach my $i (keys %rpf_total)  {
	if ( exists $rpf_cod_num{$i} )   {
		$rpf_avg_gene_dens{$i} = $rpf_total{$i}/$rpf_cod_num{$i};  }  }

###############################
#calculate coverage for each codon
##########################
my %codon_pool; # count the codon occurance after occupied codon

open (out1, ">$sample/$sample\_site_$site\_occupancy_norm1.txt");

print out1 "sample\tGENE\tPosition\tcodon\treads\tAVG_density\n";

foreach my $id (keys %si)    {

if ( exists $gene{$id} )
  {
    for (my $i=100+$ATG_distance*3+$start; $i<($si{$id}-100-$stop_distance*3-$start); $i+=3)              {
  
      
       if ( exists $rpf{$id."\t".$i} and exists $rpf_avg_gene_dens{$id} ) {
	    
       	      my $cod=substr($gene{$id},$i,3);
	      my $pool = substr($gene{$id},$i+3,9);
		
		foreach (my $j = 0; $j < 9; $j += 3)   {
			$codon_pool{substr($pool, $j, 3)} ++;  }

	      my $reads = $rpf{$id."\t".$i};

 	      my $norm = $rpf{$id."\t".$i}/$rpf_avg_gene_dens{$id};

	      print out1 "$sample\t$id\t$i\t$cod\t$reads\t$norm\n";    }   }
     }      }

close out1;

open (out2, ">$sample/$sample\_site\_$site\_occupancy_pool.txt");

foreach my $i (keys %codon_pool)    {
	print out2 "$i\t$codon_pool{$i}\n";
                                       }

close out2;    }

