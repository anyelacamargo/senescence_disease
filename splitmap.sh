#!/bin/sh


#netopt=$HOME/transsyswork/modeldiscrimination/trunk/netopt

function do_run()
{
  echo $*
  if $* ; then
    true
  else
    exit 1
  fi
}

function checkpython()
{
  if [ -z "$PYTHONPATH" ]
    then
    export PYTHONPATH=$HOME/lib64/python
  fi
}

function cutfile()
{
	chrs=("chr1A" "chr2A" "chr3A" "chr4A" "chr5A" "chr6A" "chr7A" "chr1B" "chr2B" "chr3B" "chr4B" "chr5B" "chr6B" "chr7B" "chr1D" "chr2Dx" "chr3D" "chr4D" "chr5D" "chr6D" "chr7D")
	for i in "${chrs[@]}"
	do
		echo $i
                file_name=`printf '%s.wheat.map' ${i}`
		#echo $file_name
		grep ${i} map.wheat.txt > $file_name
		echo -e "marker\tchromosome\tbp" | cat - $file_name > table_with_header
		cp table_with_header $file_name
	done

}

function dummy()
{
  dummyL=1
  for i in $dummy ; do
    dummyL=`expr $dummyL + 1 `
  done
}

function dummy2()
{
  dummyL=1
  if [[ $dummy -eq 1 ]]
  then
    dummyL=`expr $dummyL + 1 `
  fi
}

function dummy3()
{
  dummyL=1
  while test $dummy -le $dummyL ; do
    dummyL=`expr $dummyL + 1 `
  done
}

function createEmpiricalData()
{
  tp_name=`printf '%s' ${true_model} `
  do_run ./transsyswritesimsetOF -o  ${simgenex_model} -s ${rndseed} -N ${signal_to_noise} ${true_model}'_m00.tra' ${tp_name} ${tp_name}'_trace.txt'
}

function optimiseModelSynthetic()
{

  tp_name=`printf '%s' ${true_model} `
  for (( model=0; model<=${num_model}; model++ )) ; do
    model_name=`printf '%s_m%02d' ${tp_name} ${model}`
    do_run ./netopt -x ${tp_name}'_expr.txt' -o ${simgenex_model} -R ${num_optimisations} -g ${gradientfile} -T ${transformerfile} -L ${logfile} -s ${rndseed} -c ${model_name}
  done
}



function write_cshscript_header ()
{
  scriptname=$1
  echo '#!/bin/csh' > $scriptname
  echo '#$ -q long.q' >> $scriptname
  echo '#$ -j y' >> $scriptname
  echo '#$ -m e' >> $scriptname
  echo '#$ -cwd' >> $scriptname
  echo >> $scriptname
  echo 'setenv PYTHONPATH ${HOME}/lib64/python' >> $scriptname
  echo "cd $PWD" >> $scriptname
  echo 'echo start time: `date`' >> $scriptname
}


function write_cshscript_footer ()
{
  echo 'echo success: `date`' >> $scriptname
}


function add_command ()
{
  scriptname=$1
  cmd="$2"
  echo "if ( { $cmd } ) then" >> $scriptname
  echo "  echo completed \`date\`" >> $scriptname
  echo "else" >> $scriptname
  echo "  echo ERROR" >> $scriptname
  echo "  exit 1" >> $scriptname
  echo "endif" >> $scriptname
  #if [[ $rewired_topology_number -eq $num_rewired_topologies ]]
  #  then
  #  echo "qsub ${c_t}.csh" >> $scriptname
  #fi
}


function optimiseModelEmpiric()
{
  for i in $model_list ; do
    candidate_topology=${i%.*}
    scriptname=`printf 'script_%s.csh' ${candidate_topology}`
    write_cshscript_header $scriptname
    add_command $scriptname "./netopt -x  ${empirical_data} -o ${simgenex_model} -R ${num_optimisations} -g ${gradientfile} -T ${transformerfile} -L ${logfile} -s ${rndseed} -c ${candidate_topology} -r TRUE"
  do_run qsub ${scriptname}
  done
}


function mergeFile()
{
  mkdir temp
  mv opt_ja* temp
  mv tp_ja* temp
  mv *logo* temp
  cd temp
  for i in $model_list ; do
    candidate_topology=${i%.*}
    candidate_topology_logo=`printf '%s_logo.txt' ${candidate_topology}`
    cp $candidate_topology_logo $candidate_topology_logo'.bk'
    sed 's/rst/'$candidate_topology'\t/g' $candidate_topology_logo'.bk' > clean.txt
    mv clean.txt $candidate_topology_logo
    rm -f $candidate_topology_logo'.bk'
    name=`printf '%s %s ' $name $candidate_topology_logo `
  done
  fitnessname=`printf 'fitnesstable' `
  rm -fr $fitnessname'.txt'
  rm -f 'tt.txt'
  label=`printf 'model\trestart\tfitness'`
  cat $name > 'tt.txt'
  echo -e $label >> $fitnessname'.txt'
  sed '/restart/d' 'tt.txt' >> $fitnessname'.txt'
  rm -f 'tt.txt'
}


model_list=`echo *.tra`
simgenex_model=jasmonate_pipeline.sgx
signal_to_noise=0
rndseed=1
num_optimisations=6
transformerfile=transformerfile.dat
logfile=logo
gradientfile=optspec.dat
empirical_data=jasmonate_exp.txt
tp_candidate=none

while getopts m:e:s:o:d opt
do
  case "$opt" in
    m) tp_candidate="$OPTARG";;
    e) empirical_data="$OPTARG";;
    s) simgenex_model="$OPTARG";;
    o) num_optimisations="$OPTARG";;
    d) isdef=1;;
    \?) help_ani;;
  esac
done


#checkpython
#createEmpiricalData
#optimiseModelSynthetic
#optimiseModelEmpiric
#mergeFile
cutfile
