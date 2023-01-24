# first argument is genome file
# get base pair number
charcount=`grep -v ">" $1 | wc | awk '{print $3-$1}'`
echo "Base pair count is: " $charcount "(" $((charcount/1000)) " Mbp )"
totalNumberOfContigs=`grep ">" $1 | wc -l`
onePercentOfContigs=`bc -l <<< "scale=4; ($totalNumberOfContigs/100)"` 
onePercentOfContigsINT=${onePercentOfContigs%%.*}
genomeName=$1
#remove everything up to and including  last / (ie take out path)
genomeName=${genomeName##*/}
# remove .faa extention
genomeName=${genomeName::-4}
echo "Calling genome "$genomeName
#genome_chunk=$((charcount/chunkcount))
# echo "Splitting genome files into chunks of" $genome_chunk "base-pairs..."
echo "Randomly sampling genome at 5% completeness intervals pieces"
#module load bbmap
# `reformat.sh in=$1 out=./fragged_genome.faa breaklength=$genome_chunk` >> log.txt
rm -rf "fragged_"$genomeName
mkdir "fragged_"$genomeName
cd "fragged_"$genomeName
cat $1 > $genomeName".faa"
index=0
samplerate=100
while [ $index -le 18 ] # only do down to 15 % completeness in 5 % increments
do
  echo "Processing genome #" $index 
  #make new genomes from random sampling
  reformat.sh in=$genomeName".faa" out=$genomeName"_"$samplerate"_pc_complete.faa" samplerate=`bc -l <<< "$samplerate / 100"` sampleseed=1 
  # now check the actual completeness percentage
  charcount_fragged=`grep -v ">" $genomeName"_"$samplerate"_pc_complete.faa" | wc | awk '{print $3-$1}'`
  echo "Charcount_fragged: " $charcount_fragged "charcount: " $charcount 
  actualCompPercentage=`bc -l <<< "scale=6; (($charcount_fragged / $charcount) * 100)"` # we use this to label
  cat $genomeName"_"$samplerate"_pc_complete.faa" > TMP_.faa
  echo "Actual Comp is" $actualCompPercentage
  rm $genomeName"_"$samplerate"_pc_complete.faa"
  cat TMP_.faa > $genomeName"_"$actualCompPercentage"_pc_complete.faa"
  rm TMP_.faa
  

  #generate contamination
  if (($samplerate >  30)); then
    #make contamination source from incomplete genome:
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" > cont_source.faa
    f2pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 2)"`
    f5pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 5)"`
    f8pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 8)"`
    f10pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 10)"`
    f15pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 15)"`
    f20pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 20)"`
    f25pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 25)"`
    f30pc=`bc -l <<< "scale=4; ($onePercentOfContigs * 30)"`
    int2pc=${f2pc%%.*}
    int5pc=${f5pc%%.*}
    int8pc=${f8pc%%.*}
    int10pc=${f10pc%%.*}
    int15pc=${f15pc%%.*}
    int20pc=${f20pc%%.*}
    int25pc=${f25pc%%.*}
    int30pc=${f30pc%%.*}
    reformat.sh in=cont_source.faa out=2pc.faa samplereadstarget=$int2pc 
    reformat.sh in=cont_source.faa out=5pc.faa samplereadstarget=$int5pc 
    reformat.sh in=cont_source.faa out=8pc.faa samplereadstarget=$int8pc 
    reformat.sh in=cont_source.faa out=10pc.faa samplereadstarget=$int10pc 
    reformat.sh in=cont_source.faa out=15pc.faa samplereadstarget=$int15pc 
    reformat.sh in=cont_source.faa out=20pc.faa samplereadstarget=$int20pc 
    reformat.sh in=cont_source.faa out=25pc.faa samplereadstarget=$int25pc 
    reformat.sh in=cont_source.faa out=30pc.faa samplereadstarget=$int30pc 
    rm cont_source.faa

    # now calculate the size of contamination chunks
    f2pcCont=`grep -v ">" 2pc.faa | wc | awk '{print $3-$1}'`
    f2pcContActual=`bc -l <<< "scale=6; (($f2pcCont / $charcount) * 100)"` # actual basepairs of '2 percent contamination'
    
    f5pcCont=`grep -v ">" 5pc.faa | wc | awk '{print $3-$1}'`
    f5pcContActual=`bc -l <<< "scale=6; (($f5pcCont / $charcount) * 100)"` # actual basepairs of '5 percent contamination'
    
    f8pcCont=`grep -v ">" 8pc.faa | wc | awk '{print $3-$1}'`
    f8pcContActual=`bc -l <<< "scale=6; (($f8pcCont / $charcount) * 100)"` # actual basepairs of '8 percent contamination'
    
    f10pcCont=`grep -v ">" 10pc.faa | wc | awk '{print $3-$1}'`
    f10pcContActual=`bc -l <<< "scale=6; (($f10pcCont / $charcount) * 100)"` # actual basepairs of '10 percent contamination'
    
    f15pcCont=`grep -v ">" 15pc.faa | wc | awk '{print $3-$1}'`
    f15pcContActual=`bc -l <<< "scale=6; (($f15pcCont / $charcount) * 100)"` # actual basepairs of '15 percent contamination'
    
    f20pcCont=`grep -v ">" 20pc.faa | wc | awk '{print $3-$1}'`
    f20pcContActual=`bc -l <<< "scale=6; (($f20pcCont / $charcount) * 100)"` # actual basepairs of '20 percent contamination'
    
    f25pcCont=`grep -v ">" 25pc.faa | wc | awk '{print $3-$1}'`
    f25pcContActual=`bc -l <<< "scale=6; (($f25pcCont / $charcount) * 100)"` # actual basepairs of '20 percent contamination'
    
    f30pcCont=`grep -v ">" 30pc.faa | wc | awk '{print $3-$1}'`
    f30pcContActual=`bc -l <<< "scale=6; (($f30pcCont / $charcount) * 100)"` # actual basepairs of '20 percent contamination'

    
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" 2pc.faa  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f2pcContActual"pc_contaminated.faa"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" 5pc.faa  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f5pcContActual"pc_contaminated.faa"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" 8pc.faa  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f8pcContActual"pc_contaminated.faa"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" 10pc.faa  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f10pcContActual"pc_contaminated.faa"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" 15pc.faa  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f15pcContActual"pc_contaminated.faa"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" 20pc.faa  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f20pcContActual"pc_contaminated.faa"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" 25pc.faa  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f25pcContActual"pc_contaminated.faa"
    cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" 30pc.faa  > $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f30pcContActual"pc_contaminated.faa"
    rm 2pc.faa
    rm 5pc.faa
    rm 8pc.faa
    rm 10pc.faa
    rm 15pc.faa
    rm 20pc.faa
    rm 25pc.faa
    rm 30pc.faa
    rm cont_source.faa
  fi
  cat $genomeName"_"$actualCompPercentage"_pc_complete.faa" > $genomeName"_"$actualCompPercentage"_pc_completeXXX0pc_contaminated.faa"
  rm $genomeName"_"$actualCompPercentage"_pc_complete.faa"
  rm shredded_whole.faa
  if (($samplerate <  40)); then
    rm $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f30pcContActual"pc_contaminated.faa"
    rm $genomeName"_"$actualCompPercentage"_pc_completeXXX"$f25pcContActual"pc_contaminated.faa"
  fi
  index=$((index+1))
  samplerate=$((samplerate-5))
done
rm $genomeName".faa"
for i in *.faa; do sed '/^>/ s/ .*//' -i $i; done
for i in *.faa; do awk '/^>/{$0=$0"_"(++i)}1' $i > tmp"$i" && mv tmp"$i" $i; done
#mkdir ../all_genomes
#cp *.faa ../all_genomes
