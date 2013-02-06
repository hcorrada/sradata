OUTDIR="$HOME/sgseq/src/intmodel/sra"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$OUTDIR/lib64:/thumper/ctsa/genomicsR/extras/libxml2/lib"
SRA_SDK="$OUTDIR/sra_sdk-1.0.0-b10"

echo "gcc -o spot_count -I$SRA_SDK/itf spot-count.c -L$OUTDIR/lib64 -lsradb -lvdb -lklib -lkascii -lm -lz -lbz2 -ldl -lpthread"
eval "gcc -o spot_count -I$SRA_SDK/itf spot-count.c -L$OUTDIR/lib64 -lsradb -lvdb -lklib -lkascii -lm -lz -lbz2 -ldl -lpthread"

./spot_count /amber1/archive/sgseq/data/1kgenomes/low_coverage/HG00116/SRR035464
