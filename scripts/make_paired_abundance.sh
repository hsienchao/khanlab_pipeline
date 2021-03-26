abundance=$1
perl -lane '{if ($F[0] eq "pbid") {print;} else {print "$F[0].1\t$F[1]\t$F[2]\n$F[0].2\t$F[1]\t$F[2]"}}' $1 | grep -v '^#'
