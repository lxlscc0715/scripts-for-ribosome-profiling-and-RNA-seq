use strict;
use warnings;
use Number::Format qw(round);

sub process_wig  {

	my ($sam, $total) = @_;
	my $sample = $sam;
	$sample =~ s/\_\d$//;
	print " sample = $sam; $sample\n";
	my %ref_RNA;
	open (hand1, "$sample\_isoforms.fpkm_tracking") or die $!;

	while(<hand1>)  {

		chomp;
		next if /^tracking_id/;
		my @a = split /\t/;
		$ref_RNA{$a[0]} = $a[9];
		print "$a[0]\t$a[9]\n";
      }
	close hand1;	

=head1	
	my $sum;
	open (hand1, "$sam/$sam\_A.wig") or die $!;
	while (<hand1>)  {

		chomp;
		my @a = split /\t/;
		$sum += $a[3];
		}

	close hand1;

	print "total = $sum\n";
=cut

	open (hand1, "$sam/$sam\_A.bg") or die $!;
	open (out1, ">$sam/$sam\_A_normalized_total_reads.wig");
	while (<hand1>)  {

		chomp;
		my @a = split /\t/;
		if (exists $ref_RNA{$a[0]} & $ref_RNA{$a[0]} > 0 )  {
			my $value = $a[3]*100000000/($total);
			$value = round($value, 3);
			print out1 "$a[0]\t$a[1]\t$a[2]\t$value\n";
			}  }

	close hand1; close out1;
}
process_wig("RPF_sample_name1", xxxxx);    #xxxxx is total reads number
process_wig("RPF_sample_name1", xxxxx);    #xxxxx is total reads number
