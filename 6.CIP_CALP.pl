###Parsing blast results with CIP and CALP paremeters for conserved gene pairs for ancestral karyotype study###
###usage: $0 blast.output > CIP_CALP.out###
use Bio::SearchIO;

my $infile  = shift;

my $parser = Bio::SearchIO->new(-file   => $infile,
                                -format => 'blast');
while( my $result = $parser->next_result ) {
    print "Query: ",$result->query_name,"\tlength: ",$result->query_length,
          "\tnumber of hits: ",$result->num_hits,"\n";
    while(my $hit=$result->next_hit){
        print "\tHit: ",$hit->name,"\tlength: ",$hit->hit_length,
                  "\te-value: ",$hit->significance,"\n";
        my $al=0;
        my @id=();
        while(my $hsp=$hit->next_hsp){
            print "\t\tHsp length: ",$hsp->hsp_length,"\tID length: ",
                          $hsp->num_identical,"\t idenity: ",
                          $hsp->percent_identity,"\n";
            $al=$al+($hsp->hsp_length);
            my $ai=($hsp->num_identical)/($hsp->hsp_length);
            push @id,$ai;
        }
        my $cip=0;
        foreach (@id){
            $cip=$cip+($_/$al*100);
            }
        my $calp=$al/($result->query_length);
        print "\tcip: ",$cip,"\tcalp: ",$calp,"\n";
    }
}
