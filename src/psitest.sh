#!/bin/bash

. ~/src/common/util.sh

usage () {
    echo "$(basename $0) -g|p"
    exit 2
}

niters=1001 ## was 1000 but termination condition in em.c had off-by-one bug
nprocs=1
while getopts "gpdP:" opt ; do 
    case $opt in
	p) dataarg="-p $HOME/src/structure/test/genotypes.chiamo" ;;
	g) dataarg="-g $HOME/src/structure/test/genotypes.geno" ;;
	d) dummy=$true ;;
	P) nprocs=$OPTARG ;;
	*) usage ;;
    esac
done

[ -n "$dataarg" ] || usage

[ -n "$dummy" ] && {
    outdir=tmp.out 
    rm -rf $outdir
} || outdir=~/src/structure/test/out/$(date +%Y-%m-%d-%s)

mkdir $outdir
tmpfile=/var/tmp/psi-test

if [ $HOSTNAME = zuse.osc.ox.ac.uk ] ; then
    export NSLOTS
    export TMPDIR
fi

/usr/bin/time -f %e -o $tmpfile \
    ~/src/structure/psi.sh \
    -n 80 -L 1000 -K 2 -t $niters \
    $dataarg \
    -m ~/src/structure/test/mu-init \
    -q ~/src/structure/test/q-init \
    -P $nprocs \
    -o $outdir \


[ $? ] || {
    echo error!
    rm -r $outdir
    exit 2
}
echo -n "$(basename $outdir) $dataarg   " | cat - $tmpfile > $outdir/time

echo dataarg: $dataarg
echo log-posteriors:
cat $outdir/loglike-$(printf '%05d' $niters)
echo -n "should be" ; cat ~/src/structure/test/logpost

echo times:
head -q -n1 ~/src/structure/test/out/*/time

exit 0
