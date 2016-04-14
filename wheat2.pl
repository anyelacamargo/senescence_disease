#!/usr/bin/perl
use strict;

my $genofile = "wheat_geno2.csv";
open(GENOS,$genofile) || die "Could not open $genofile\n";

my $file = "founders.geno2.txt";
open(FILE,$file) || die "Could not open $file\n";

my $mapfile = "map.wheat.txt";
open(MAPFILE, ">$mapfile") || die "Could not open $mapfile\n";
print MAPFILE "marker	chromosome	bp\n";

$_ = <FILE>;
chomp;
my ( $chrno, $chr, $pos, $posjittered, $genomewidepos, $name, $X13074, $X , @founders ) = split;
my $geno={};
my $data={};
my $map={};
my $use={};
my $fped = {};

my $nf = $#founders+1;
print "@founders\n";
while(<FILE>) {
  chomp;
  my ( $N, $chrno, $chr, $pos, $posjittered, $genomewidepos, $name, $X13074, @genos ) = split;
  $data->{$X13074} = \@genos;
  push @{$map->{$chr}}, [$chr, $posjittered, $X13074];
}

my @na = ();
for( my $k=0;$k<$nf;$k++) {
  $na[$k] = sprintf( "%.4f", 1/$nf );
}

foreach my $chr ( sort keys %$map ) {
  my $bp=1;
  my $allelesfile = "chr$chr.wheat.alleles";
  my $allelestxt = "";
  my $pedtxt = "";

  foreach my $f ( @founders ) {
    $fped->{$chr}{$f} = "founder\t$f\t0\t0\t1\t0";
  }
  my @mark = @{$map->{$chr}};
  my $nm = 0;
  foreach my $ref ( @mark ) {
    my ( $chr, $pos, $snp ) = @$ref;
    my @g = @{$data->{$snp}};
    my ( @a0, @a1, $n0, $n1);
    my @b;
    for( my $k=0;$k<=$#g;$k++) {
      if ( $g[$k]==0 ) {
	push @a0, 1;
	push @a1, 0;
	++$n0;
	$b[$k] = 0;

      } else {
	push @a1, 1;
	push @a0, 0;
	++$n1;
	$b[$k] = 1;
      }
    }
    if ( $n0*$n1 > 0 ) {
      $allelestxt .= "marker $snp 3 $chr $pos\n";
      $nm++;
      $use->{$snp} = 1;
      for( my $k=0;$k<=$#g;$k++) {
	$fped->{$chr}{$founders[$k]} .= "\t$b[$k] $b[$k]";
	$a0[$k] = sprintf( "%.4f", $a0[$k]/$n0 );
	$a1[$k] = sprintf( "%.4f", $a1[$k]/$n1 );
      }
      
      $allelestxt .= "allele 0 @a0\n";
      $allelestxt .= "allele 1 @a1\n";
      $allelestxt .= "allele NA @na\n";
      $bp += 100000;
      print MAPFILE "$snp\tchr$chr\t$bp\n";
    }
  }

  open(ALLELES, ">$allelesfile" ) || die "Could not open $allelesfile\n";
  print ALLELES "markers $nm strains $nf\nstrain_names @founders\n$allelestxt";
  print "writing $allelesfile\n";
  close(ALLELES);
}
my $h1 = <GENOS>;

chomp;
$h1 = <GENOS>;

chomp;
my ( $a, $b, @h1) = split(/,/, $h1 );
$_ = <GENOS>;
s/\r\n|\n|\r/\n/g;
chomp;
my ( $a, $b, @h2) = split(/,/);

while(<GENOS>) {
  s/AA/0 0/g;
  s/BB/1 1/g;
  s/AB/0 1/g;
  s/NC/NA NA/g;
  s/Null/NA NA/g;
  s/\r\n|\n|\r/\n/g;
  chomp;
 
  my ( $a, $snp, @g )  = split( /,/ );
  if ( $use->{$snp} ) {
    $geno->{$snp} = \@g;
  }
}
close(GENOS);
foreach my $chr ( sort keys %$map ) {
  my $pedfile = "chr$chr.wheat.data";
  open(PED, ">$pedfile") || die "Could not open $pedfile\n";
  print "writing $pedfile\n";
  my $data = "";
  my @cmap = @{$map->{$chr}};
  my @snps;
  foreach my $ref ( @cmap ) {
    my $snp = $ref->[2];
    push @snps, $snp if ( $use->{$snp} );
  }
  
  for( my $k=0;$k<=$#h1;$k++) {
    $h1[$k] =~ s/\s+//g;
    $data .= "wheat\t$h1[$k]\t0\t0\t1\t0\t";
    my @g = ();
    foreach my $snp ( @snps ) {
      my $g = $geno->{$snp}[$k] ? $geno->{$snp}[$k] :"NA NA";
      push @g, $g;
    }
    $data .= join( "\t", @g ) . "\n";
  }
  foreach my $f ( @founders ) {
    print PED $fped->{$chr}{$f} . "\n";
  }
  print PED $data;
  close(PED);
}

close(MAPFILE);


