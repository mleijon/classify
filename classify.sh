#!/bin/bash

shopt -s nocaseglob
export PATH=$PATH:$HOME/classify


#######INPUT DATA CONTROL########
max_size=${2:-50000000}
if [ ! -d "$1" ]; then
  echo "Enter the full path to the directory containing the sequencing \
data as argument."
  exit 0
else
  FILES=($1/*R1*.gz)
  if [ ! -f ${FILES[0]} ]; then
    echo "No sequence files in: $1"
    exit
  fi
  dir=$(dirname "$FILES")
fi
if [[ "${2:0:-1}" =~ ^[0-9]+$ ]]; then
  :
elif [ -z $2 ]; then
  :
else
  echo "Input error: $2"
  exit
fi
if [ ! -z $2 ]; then
  case ${2: -1} in
    G|g)
      max_size=$((1000000000*${2:0:-1}))
      ;;
    M|m)
      max_size=$((1000000*${2:0:-1}))
      ;;
    [0-9])
      ;;
    *)
      echo "Input error: $2"
      exit
      ;;
  esac
fi
##################################
OUTDIR=$dir/classification
if [ -d $OUTDIR ]; then
  rm -r $OUTDIR
fi
mkdir $OUTDIR
for f in ${FILES[@]}; do
  base=$(basename "$f")
  ext=${f##*.}
      if [ -d $OUTDIR/${base%%_*} ]; then
    rm -r $OUTDIR/${base%%_*}
  fi
  mkdir $OUTDIR/${base%%_*}
  nr_lines=$(gunzip -c $f|wc -l)
  nr_reads=$((nr_lines/4))
  file_size_fwd=$(wc -c $f)
  file_size_fwd=${file_size_fwd% *}
  FILE=${f/R1/R*}
  FILE=${FILE##*/}
  if [ ! -f ${f/R1/R2} ]; then
    echo "No reverse file: ${f/R1/R2}. Skipping this sample."
    continue
  fi
  file_size_rev=$(wc -c ${f/R1/R2})
  file_size_rev=${file_size_rev% *}
  if [ $file_size_fwd -lt $file_size_rev ]; then
    file_size=$file_size_rev
  else
    file_size=$file_size_fwd
  fi

  if [ $file_size -le $max_size ]; then
    echo "No splitting required for: ${f/R1/R*}. Copying files..."
    if [ -d $OUTDIR/${base%%_*} ]; then
      rm -r $OUTDIR/${base%%_*}
    fi
    mkdir $OUTDIR/${base%%_*}
    cp $f $OUTDIR/${base%%_*}
    cp ${f/R1/R2} $OUTDIR/${base%%_*}
  else
    if [ ! $((file_size%max_size)) -eq 0 ]; then
      nr_files=$((file_size/max_size + 1))
    else
      nr_files=$((file_size/max_size))
    fi
    if [ $nr_files -gt 100 ]; then
      echo "Too many splits: $nr_files (maximum: 100)"
      exit
    fi
    split_size=$((4*(nr_reads/nr_files + nr_reads%nr_files)))
    as='_'${base#*_}
    echo "Splitting $FILE..."
    gunzip -c $f|split -l $split_size --filter='gzip > $FILE.gz' \
    --additional-suffix=${as%.$ext} - ${f%%_*}'xxx'
    rev=${f/R1/R2}
    as=${as/R1/R2}
    gunzip -c ${f/R1/R2}|split -l $split_size --filter='gzip > $FILE.gz' \
     --additional-suffix=${as%.$ext} - ${rev%%_*}'xxx'
    mv $dir/*xxx* $OUTDIR/${base%%_*}
  fi

  cp $HOME/PycharmProjects/classify/daa2spec.py $HOME/classify/
  cp $HOME/PycharmProjects/classify/rundiam_lf.sh $HOME/classify/
  #RUNDIAM_LF
  echo "Processing $FILE..."
  wait;rundiam_lf.sh $OUTDIR/${base%%_*}
  rm $OUTDIR/${base%%_*}/*.gz
  #cat $OUTDIR/${base%%_*}/*/*.gz >> $OUTDIR/${base%%_*}_uq.fasta.gz
  cat $OUTDIR/${base%%_*}/*/*.daa >> $OUTDIR/${base%%_*}.daa
  #rm -r $OUTDIR/${base%%_*}

  #DAA2SPEC.PY
  echo "Classifying [daa2spec.py]"
  wait;daa2spec.py -f $OUTDIR/${base%%_*}.daa -a -b -s -e -u -o -d -v --derep\
  &>/dev/null
  gzip $OUTDIR/${base%%_*}/*/*.daa
  gzip $OUTDIR/${base%%_*}.daa
done



