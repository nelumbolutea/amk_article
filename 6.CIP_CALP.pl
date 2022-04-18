#####filter blastp output for ancestral karyotype analysis by CIP & CALP parameters###
#####Both CIP CALP should >=0.5###
##e.g.##
##makeblastdb -in species2.pep  -dbtype prot -out species2.pep.db##
##blastp -query species1.pep -db species2.pep.db -num_threads 12 -evalue 1E-5 -outfmt 0  -out species1-species2.blastout###
##perl cip_calp.pl species1-species2.blastout>species1-species2.blastout.filtered###

use Bio::SearchIO;

my $blastout  = shift;

my $filter = Bio::SearchIO->new(-file   => $blastout,
                                -format => 'blast');
print "Query\tlength_Q\tnumber_of_hits\tHit\tlength_H\te-value\tCIP\tCALP\n";				  
while( my $result = $filter->next_result ) {
    
    while(my $hit=$result->next_hit){
	 
        my $al=0;
        my @id=();
        while(my $hsp=$hit->next_hsp){

            $al=$al+($hsp->hsp_length);
            my $ai=($hsp->num_identical)/($hsp->hsp_length);
            push @id,$ai;
        }
        my $cip=0;
        foreach (@id){
            $cip=$cip+($_/$al*100);
            }
        my $calp=$al/($result->query_length);
	  if($cip ge 0.5 and $calp ge 0.5)
	  {
	  print $result->query_name,"\t",$result->query_length,"\t",$result->num_hits,"\t",$hit->name,"\t",$hit->hit_length,"\t",$hit->significance,"\t",$cip,"\t",$calp,"\n";	
	  }
    }
}
