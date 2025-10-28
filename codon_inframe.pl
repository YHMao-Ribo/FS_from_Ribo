use strict;use warnings;

my $cdna=shift;
my $reads=shift;

my %cdna;
open(F,"gzip -dc $cdna|");
while(my $line=<F>){
	$line=~s/\s+$//;$line=~s/>//;
	my $seq=<F>;$seq=~s/\s+$//;

	$cdna{$line}=$seq;
}
close(F);

open(F,"gzip -dc $reads|");
my @reads=<F>;
close(F);


my @shuffle_reads;
for(my $i=0;$i<scalar(@reads);$i++){
	my $rand_idx=int(rand(scalar(@reads)-$i))+$i;
	($reads[$i],$reads[$rand_idx])=($reads[$rand_idx],$reads[$i]);
	push @shuffle_reads,$reads[$i];
}
my $break=int(scalar(@shuffle_reads)/5);

my %mean;my %mean_count;
my %frame_mRNA;
for(my $i=0;$i<scalar(@shuffle_reads);$i++){
	my $line=$shuffle_reads[$i];
	$line=~s/\s+$//;
	my @tps=split(/,/,$line);
	my $name=shift @tps;
	
	my($start,$stop,$tr_len)=$name=~/(\d+):(\d+):(\d+$)/;

	for(my $i=$start+150;$i<$stop-30;$i+=3){
		my $codon=substr($cdna{$name},$i,6);

		my $stop_pos_frame2=-1;
		for(my $j=$i+2;$j<$tr_len-3;$j+=3){
			my $stop=substr($cdna{$name},$j,3);
			if($stop eq "TAA" or $stop eq "TGA" or $stop eq "TAG"){
				$stop_pos_frame2=$j;last;
			}
		}

		#		$stop_pos_frame2=$i+60+2 if($stop_pos_frame2-$i>60);

		if($stop_pos_frame2>0){
			for(my $j=$i-30;$j<$stop_pos_frame2;$j+=3){
				next if($i<0 or $j>$stop-3);
				next if($tps[$j]+$tps[$j+1]+$tps[$j+2]>1000);

				if($j<$i){
				}
				else{
					$frame_mRNA{$codon}{'inframe'}+=$tps[$j];
					$frame_mRNA{$codon}{'total'}+=($tps[$j]+$tps[$j+1]+$tps[$j+2]);
				}	
			}
		}
	}

	if($i>0 and $i%($break+1)==0){
		foreach my $codon(keys %frame_mRNA){
			if(exists $frame_mRNA{$codon}{'total'} and $frame_mRNA{$codon}{'total'}>5){
				$frame_mRNA{$codon}{'inframe'}/=$frame_mRNA{$codon}{'total'};
				$mean{$codon}{'inframe'}+=$frame_mRNA{$codon}{'inframe'};
				$mean{$codon}{'total'}++;
			}
		}
		%frame_mRNA=();
	}
}
close(F);

foreach my $codon(keys %frame_mRNA){
	if(exists $frame_mRNA{$codon}{'total'} and $frame_mRNA{$codon}{'total'}>5){
		$frame_mRNA{$codon}{'inframe'}/=$frame_mRNA{$codon}{'total'};
		$mean{$codon}{'inframe'}+=$frame_mRNA{$codon}{'inframe'};
		$mean{$codon}{'total'}++;
	}
}

$reads=~s/^\S+\///;
open(F,">2double_codon_IFR_$reads");
foreach my $codon(sort keys %mean){
	my($codon1,$codon2)=$codon=~/(\w{3})(\w{3})/;
	$codon1=join(",",$codon1,$codon2);
	next if($codon1=~/TAA/ or $codon1=~/TAG/ or $codon1=~/TGA/);
	
	$mean{$codon}{'inframe'}/=$mean{$codon}{'total'};
	print F "$codon\t$mean{$codon}{'inframe'}\t$mean{$codon}{'total'}\n";
}
close(F);
