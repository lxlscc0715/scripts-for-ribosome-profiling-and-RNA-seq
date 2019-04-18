use strict;
use warnings;

=head1
This script is used to remove adaptor sequences, specifically for ribosome profiling or sRNA seq, etc. Please carefully check the adaptor sequences before use it.
=cut

my $min=12;
my $max=150;

&trim3 ('sample_name1');
&trim3 ('sample_name2');
&trim3 ('sample_name3');
&trim3 ('sample_name4');

####trim 3'
sub trim3                                                      {
my ($sam)=@_; my $count=0; my $id=10000000; my $idfq; my $on=0; my $si=0;
my $cca=0;my %cca=();

mkdir $sam;
mkdir 'fq';
open (hand1,"fq/$sam.fq") or die $!;
open (hand2,">$sam/$sam\_$min.txt");
open (hand3,">fq/$sam\_$min.fq");

while (<hand1>)                 {
$_ =~ s/\s+$//;
$count++;
my $mod=$count % 4;

if ($mod == 1)       {
$idfq=$_; next  }

if ($mod == 2)             {
$on=0; $si=0;
if (/^(.*[ATCG]{4}AGATCGGAAG)/)   {
$_=$1;
$_ =~ s/[ATCG]{4}AGATCGGAAG$//;   }
else                 {
$_ =~ s/([ATCG]{4}AGATCGGA|[ATCG]{4}AGATCGG|[ATCG]{4}AGATCG|[ATCG]{4}AGATC|[ATCG]{4}AGAT)[ATCGN]*$//;
                     }

if (length($_)>$min and length($_)<$max)   {
$on=1; $id++;
$si=length($_);
print hand2 $_,"\n";
print hand3 "$idfq\n";
print hand3 $_,"\n";
print hand3 "+$id\n";}     }
elsif ($mod ==0 && $on==1) {
my $qua=substr ($_,0,$si);
print hand3 $qua,"\n";     }   }
close hand1; close hand2; close hand3;

open (hand1,"$sam/$sam\_$min.txt");
my %uni=(); $count=0; my %sta=(); my %siz=(); my $number=0;
while (<hand1>)            {
$_ =~ s/\s+$//; $count++;
$uni{$_}++;
$sta{substr($_,0,1)}++;
$siz{length($_)}++;        }
$number=scalar keys %uni;
close hand1;
unlink "$sam/$sam\_$min.txt";

$id=10000000;
open (hand1,">$sam/$sam\_$min\_uni.txt");
foreach my $seq (sort {$uni{$b}<=>$uni{$a}} keys %uni)  {
$id++;
print hand1 ">$id\_x$uni{$seq}\n$seq\n";                }
close hand1; %uni=();

open (hand1,">$sam/sum\_$sam.txt");
print hand1 "$sam\n\nTotal number of $min nt or longer reads\t$count\n";
print hand1 "Total unique species of $min nt or longer\t$number\n";
close hand1;

open (hand1,">>$sam/sum\_$sam.txt");
print hand1 "\nstart:\n";
foreach my $st (keys %sta)           {
print hand1 $st,"\t",$sta{$st},"\n"; }
%sta=();

print hand1 "\nsize:\n";
foreach my $si (sort {$a<=>$b} keys %siz)   {
print hand1 $si,"\t",$siz{$si},"\n";        }

close hand1; %siz=(); %sta=();                                   }
