package MyBioinfo::NBTest;

use strict;
use Math::CDF;
use MyBioinfo::Common;

require Exporter;
our @ISA = qw(Exporter);

our @EXPORT_OK = qw();

our @EXPORT = qw(nb_pval);

our $VERSION = '0.50';

#r,p parameters of negative binomial distribution
sub nb_r_p{
	my ($m,$v,$ref,$eps)=@_;
	my $r;
	my $p;
#	print "v = $v\t m = $m\n";
	if($m==0)
	{
		$r = $eps;
		$p = $eps;
	}
	else{
	 
	 if($v<$m)
	 {
	 	$v = $m + $eps;
	 }	
	 $r = $m**2/($v-$m);
	 $p = $r/($r + $m);
	}
	push @{$ref},$r;
	push @{$ref},$p;		
}

#density distribution of negative binomial 
sub dnbinom{
	my ($x,$r,$p)=@_;
	#print "$x\t$r\t$p\n";
	
	if($x-1<0)
	{
		if($x==0)
		{
			return &Math::CDF::pnbinom($x,$r,$p) +1e-08;			
		}
		else 
		{		
			return 	&Math::CDF::pnbinom($x,$r,$p)-&Math::CDF::pnbinom(0,$r,$p) +1e-08;
			
		}
		
	} 
	else{
		return  &Math::CDF::pnbinom($x,$r,$p)-&Math::CDF::pnbinom($x-1,$r,$p) + 1e-08; #avoid zero case	
	}	
	
}

#quick computation of pvalues of negative binomial distribution
sub nb_join_pval{
	my($kl,$kr,$kS,$sizeA,$probA,$sizeB,$probB,$probs,$eps) = @_;
	 
    my $lval = dnbinom( $kl, $sizeA, $probA ) * dnbinom( $kS-$kl, $sizeB, $probB );
    my $rval = dnbinom( $kr, $sizeA, $probA ) * dnbinom( $kS-$kr, $sizeB, $probB );
    my $prevlval = $lval;
    my $prevrval = $rval;
    my $total = $lval + $rval;
    my $esttotalperlength = $total/2;
    
    my $obstotal = 0;   
    my $step = 1;
    my $steps = 0;
    my $do_left ;
    
   if( $lval <= $probs )
   {   $obstotal += $lval; }
   if( $rval <= $probs )
   {   $obstotal += $rval; }
   
   while( $kl < $kr ) {
      $steps ++;
      if( abs($prevrval - $rval) / $prevrval > .01 )
      {   $do_left = 1;}
      else{ if( abs($prevlval - $lval) / $prevlval > .01 )
         { $do_left = 0; }
      else
         {$do_left = $lval > $rval;}
      }
            
      if( $do_left ) { #left to right
         $prevlval = $lval;
         if( $kl + $step > $kr )
         {   $step = $kr - $kl;}
         $kl += $step;
         $lval = dnbinom( $kl, $sizeA, $probA) * dnbinom( $kS-$kl, $sizeB, $probB);
      
         if( $step == 1 )
         {   $total += $lval; }   
         else
         {   $total += min( $lval, $prevlval ) * $step; }   #estimate multiple steps
         
         if( $lval <= $probs ) {
            if( $step == 1 )
            {   $obstotal += $lval;}
            else {       
               if( $prevlval <= $probs )
               {   $obstotal += max( $lval, $prevlval ) * $step;}
               else
               {   $obstotal += max( $lval, $prevlval ) * $step * abs( ($probs-$lval) / ($prevlval-$lval) );}
            }
         }       
         if( abs( $prevlval - $lval ) / $prevlval < $eps )
         {   $step = max( $step + 1, $step * 1.5 );}
            
      } 
      
      else { #right to left
         $prevrval = $rval;
         if( $kr - $step < $kl )
          {  $step = $kr - $kl;}
         $kr -= $step;
         $rval = dnbinom( $kr, $sizeA, $probA ) * dnbinom( $kS-$kr, $sizeB, $probB);
         if( $step == 1 )
         {   $total += $rval;}   
         else
         {   $total += min($rval, $prevrval ) * $step;}   
         if( $rval <= $probs ) {
            if( $step == 1 )
             {  $obstotal += $rval;}
            else {       
               if( $prevrval <= $probs )
               {   $obstotal += max( $rval, $prevrval ) * $step;}
               else
               {   $obstotal += max( $rval, $prevrval ) * $step * abs($probs-$rval) / ($prevrval-$rval);}
            }
         }       
         if( abs( $prevrval - $rval) / $prevrval < $eps )
         {   $step = max( $step + 1, $step * 1.5 );}
      }
   }
   
   if($obstotal > $total)
   {
   	return ($total,$total);
   }
  return ($total, $obstotal);
} 

#p-value of negative binomial distribution
sub nb_pval{
	my ($ka,$kb,$miu1,$var1,$miu2,$var2,$eps) = @_;
	
	my $ks = $ka+$kb;
	
	my @rp1 = ();
	nb_r_p($miu1,$var1,\@rp1,$eps);
	my @rp2 = ();
	nb_r_p($miu2,$var2,\@rp2,$eps);
	
	my $r1 = $rp1[0];
	my $p1 = $rp1[1];
	
	my $r2 = $rp2[0];
	my $p2 = $rp2[1];
	
	
#	print "$miu1,$var1\n";
	#print "A = $ka\t$r1\t$p1\n";
	#print "B = $kb\t$r2\t$p2\n";
	my $kexp = ($ka+$kb) * $miu1/($miu1 + $miu2); 
	
	my $pv1 = dnbinom($ka,$r1,$p1);
	my $pv2 = dnbinom($kb,$r2,$p2);
	my $p = $pv1 * $pv2;
	#print "$pv1\t$pv2\t$p\n";
	
	my ($ldenom,$lnume) = nb_join_pval(0,$kexp,$ka+$kb,$r1,$p1,$r2,$p2,$p,$eps);
	
	my $kl2 = $kexp +1;
	if($kexp +1> $ka+$kb)
	{
		$kl2 = $ka + $kb;
	}
	my ($denominator,$numerator) = nb_join_pval($kl2,$ka+$kb,$ka+$kb,$r1,$p1,$r2,$p2,$p,$eps);
	if($numerator + $lnume==0)
	{
		
		return $p/($denominator+$ldenom);
	}
	return ($numerator+$lnume)/($denominator+$ldenom);	
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

MyBioinfo::NBTest - methods to perform negative binomial tests.

=head1 SYNOPSIS

  use MyBioinfo::NBTest;

  my $pval = &nb_pval($ka, $kb, $mu1, $var1, $mu2, $var2, $eps);  # NB-test p-value.


=head1 DESCRIPTION

  Perform NB test between two groups of counts.

=head2 EXPORT

  nb_pval

=head1 SEE ALSO

MyBioinfo::Common
MyBioinfo::Math

=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html


=cut
