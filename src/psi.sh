#!/bin/bash

. ~/src/common/util.sh

usage () {
    echo $(basename $0) "-[g|p] datafile -K K -n n -L L -o outdir -t numiters"
    exit 1
}

resume=$false

while getopts "g:p:K:n:L:o:t:m:q:P:r" opt ; do 
    case $opt in
	g) dataarg="-g $OPTARG" ;;
	p) dataarg="-p $OPTARG" ;;
	K) K=$OPTARG ;;
	n) n=$OPTARG ;;
	L) L=$OPTARG ;;
	o) outdir=$OPTARG ;;
	t) niters=$OPTARG ;;
	r) resume=$true ;;
	m) extraargs="$extraargs -m $OPTARG" ;;
	P) nprocs=$OPTARG ;;
	q) extraargs="$extraargs -q $OPTARG" ;;
	v) extraargs="$extraargs -v $OPTARG" ;;
	*) usage ;;
    esac
done
shift $(($OPTIND - 1))

[[ -n "$K" && -n "$n" && -n "$L" && -n "$outdir" && -n "$niters" && -n "$dataarg" ]] || usage

[ $resume = $true ] && {
    die "resume requested but haven't really checked to see code in psi.sh is sensible" 3
    [ -d $outdir ] || die "resume requested but outdir $outdir is not a directory" 3
    [ -n "$extraargs" ] && die "resume requested but seems like mu/q init files have been specifed separately" 3
    parentdir=$outdir
    outdir=$outdir/resume
    mkdir -p $outdir
    [ -e $outdir/mu-init ] || cat $parentdir/mu-* > $outdir/mu-init
    [ -e $outdir/q-init ] || ln -s $parentdir/q $outdir/q-init
    extraargs="-m $outdir/mu-init -q $outdir/q-init"
} || mkdir -p $outdir

logfile=$outdir/log

exec 3>&1 ## link file descriptor 3 to stdout (advanced bash guide example 19-2)
exec > $logfile ## subsequent output goes to $logfile

## work out how many processes have been requested
if [ -n "$NPROCS" ] ; then
    [ -z "$nprocs" ] || echo "warning: nprocs=$nprocs, but using NPROCS=$NPROCS inherited from qsub"
else
    [ -n "$nprocs" ] || die "error: must specify nprocs (NPROCS variable doesn't exist so presumably you are not using qsub)"
fi
nprocs=${NPROCS:-$nprocs}

echo psi: parallel structure inference
echo $(date)
echo
echo at shell level:
echo dataarg: $dataarg
echo n=$n
echo K=$K
echo L=$L
echo nprocs=$nprocs
echo MACHINEFILE=$MACHINEFILE
echo outdir=$outdir
[ $resume = $true ] && echo resumption of output in $parentdir

args="-n $n -L $L -K $K -t $niters $dataarg -o $outdir $extraargs"
echo args are $args

## command="mpirun -np $nprocs"
[ -n "$MACHINEFILE" ] && command="$command -machinefile $MACHINEFILE"
echo "going to C..." ; echo

exec 1>&3 3>&- ## restore stdout and close fd 3

$command psi $args | tee -a $logfile

exit 0
