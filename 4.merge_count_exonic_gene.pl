#! /usr/local/perl -w

open(INA,"count/Abdel2023_Upenn_NG_Ctrl1_exonic_gene_Jerry.txt")||die("Can't open INA:$!\n");
open(INB,"count/Abdel2023_Upenn_NG_Ctrl2_exonic_gene_Jerry.txt")||die("Can't open INA:$!\n");
open(INC,"count/Abdel2023_Upenn_NG_Ctrl3_exonic_gene_Jerry.txt")||die("Can't open INA:$!\n");

open(INE,"count/Abdel2023_Upenn_NG_Treat1_exonic_gene_Jerry.txt")||die("Can't open INA:$!\n");
open(INF,"count/Abdel2023_Upenn_NG_Treat2_exonic_gene_Jerry.txt")||die("Can't open INA:$!\n");
open(ING,"count/Abdel2023_Upenn_NG_Treat3_exonic_gene_Jerry.txt")||die("Can't open INA:$!\n");

open(OUT,">Merged_Abdel2023_NG_Ctrl_Treat_exonic_gene_Jerry.xls")||die("Can't write OUT:$!\n");

while(<INA>){chomp;if(/^(ENSMUSG\d+)\_[\w\.\-\(\)]+\t(\d+)$/){$ctrl1{$1}=$2;}else{print"error1\t$_\n";}}
while(<INB>){chomp;if(/^(ENSMUSG\d+)\_[\w\.\-\(\)]+\t(\d+)$/){$ctrl2{$1}=$2;}else{print"error2\t$_\n";}}
while(<INC>){chomp;if(/^(ENSMUSG\d+)\_[\w\.\-\(\)]+\t(\d+)$/){$ctrl3{$1}=$2;}else{print"error3\t$_\n";}}

while(<INE>){chomp;if(/^(ENSMUSG\d+)\_[\w\.\-\(\)]+\t(\d+)$/){$treat1{$1}=$2;}else{print"error7\t$_\n";}}
while(<INF>){chomp;if(/^(ENSMUSG\d+)\_[\w\.\-\(\)]+\t(\d+)$/){$treat2{$1}=$2;}else{print"error8\t$_\n";}}
while(<ING>){chomp;if(/^(ENSMUSG\d+)\_[\w\.\-\(\)]+\t(\d+)$/){$treat3{$1}=$2;}else{print"error9\t$_\n";}}

print OUT "\tCtrl_rep1\tCtrl_rep2\tCtrl_rep3\t";
print OUT "Treat_rep1\tTreat_rep2\tTreat_rep3\n";

my $a1=0;
foreach (keys %ctrl1)
{
 $a1++;
 print OUT "$_\t";
 print OUT "$ctrl1{$_}\t$ctrl2{$_}\t$ctrl3{$_}\t";
 print OUT "$treat1{$_}\t$treat2{$_}\t$treat3{$_}\n";
}

print "total-genes: $a1\n";

close OUT;


