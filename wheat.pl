#!/usr/bin/perl
use strict;

my $file = "wheat.txt";
open(FILE,$file) || die "Could not open $file\n";
$_ = <FILE>;
chomp;
my ( $chrno, $chr, $pos, $posjittered, $genomewidepos, $name, $X13074, $X , @founders ) = split;
my $data={};
my $map={};
my $nf = $#founders+1;
print "@founders\n";
while(<FILE>) {
  chomp;
  my ( $N, $chrno, $chr, $pos, $posjittered, $genomewidepos, $name, $X13074, $X , @genos ) = split;
  $data->{$X13074} = \@genos;
  push @{$map->{$chr}}, [$chr, $posjittered, $X13074];
}

my @na = ();
for( my $k=0;$k<$nf;$k++) {
  $na[$k] = sprintf( "%.4f", 1/$nf );
}

foreach my $chr ( sort keys %$map ) {
  my $allelesfile = "chr$chr.alleles.txt";
  my @mark = @{$map->{$chr}};
  my $nm = $#mark+1;
  open(ALLELES, ">$allelesfile" ) || die "Could not open $allelesfile\n";
  print ALLELES "markers $nm strains $nf\nstrain_names @founders\n";
  print "writing $allelesfile\n";
  foreach my $ref ( @mark ) {
    my ( $chr, $pos, $snp ) = @$ref;
    print "@$ref\n";
    my @g = @{$data->{$snp}};
    print ALLELES "marker $snp 3 $chr $pos\n";
    my ( @a0, @a1, $n0, $n1);
    for( my $k=0;$k<=$#g;$k++) {
      if ( $g[$k]==0 ) {
	push @a0, 1;
	push @a1, 0;
	$n0++;
      } else {
	push @a1, 1;
	push @a0, 0;
	$n1++;
      }      
    }
    for( my $k=0;$k<=$#g;$k++) {
      $a0[$k] = sprintf( "%.4f", $a0[$k]/$n0 );
      $a1[$k] = sprintf( "%.4f", $a1[$k]/$n1 );
    }
    
    print ALLELES "allele 0 @a0\n";
    print ALLELES "allele 1 @a1\n";
    print ALLELES "allele NA @na\n";
  }
  close(ALLELES);
}

