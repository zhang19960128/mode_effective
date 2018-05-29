#!/bin/bash
function grepnorm(){
    count=$1;
    file_name=$2;
    mass=$3;
    flag="vibration";
    atomall=${#mass[@]};
    declare -a newmode;
    grep -m$count -A $atomall "$flag" $file_name | tac | grep -m1 -B $atomall "$flag" | tac | grep -v $flag >temp1.txt;
    ####################### transform into eigenvector of dynamic matrix ###########################
    declare -a allvector
    for i in `seq 1 $atomall`
    do
        mode=`sed -n "${i}p" temp1.txt`;
        tick=`expr $i - 1`;
        mode=($mode);
        for j in `seq 0 2`
        do
            m=${mass[$tick]};
            newmode[$j]=`bc -l <<< "${mode[$j]}*sqrt($m)"`;
            allvector+=(${newmode[$j]})
        done
    done
    sum=0.0;
    for i in ${allvector[@]}
    do
        sum=`bc -l <<< "$sum+$i*$i"`;
    done
    factor=`bc -l <<< "sqrt($sum)"`
    len=${#allvector[@]};
    len=`expr $len - 1`;
    for i in `seq 0 $len`
    do
        allvector[$i]=`bc -l <<<"${allvector[$i]}/$factor"`;
    done
    ####################### end transform into eigenvector of dynamic matrix ########################
    echo ${allvector[@]}
}
function effectcharge(){
	flag_all="Effective Charges E-U: Z_{alpha}{s,beta}";
	flag="atom #";
	atom_num=$1;
	atom_all_num=$2;
	lines=$( expr $atom_all_num \* 4 );
	lines=$( expr $lines + 2 );
	file_name=$3;
	match_num="";
	for i in `seq 1 ${#atom_num}`
	do
		start=$( expr $i - 1 );
		match_num=$match_num"[""`echo ${atom_num:$start:1}`""]{1}";
	done
  grep -A $lines "$flag_all" $file_name | grep -A 3 -E "[a]{1}[t]{1}[o]{1}[m]{1}[ ]{1}[\#]{1,10}[ ]{0,7}""$match_num" | grep -v "atom"
}
function matrixtime(){
   matrix=$1;
   vector=$2;
   matrix=($matrix);
   vector=($vector);
   declare -a re=( 0 0 0 );
   for i in `seq 0 2`
   do
      for j in `seq 0 2`
      do
         tick=`expr 3 \* $i + $j`;
         re[$i]=`bc -l <<< "${re[$i]}+${matrix[$tick]}*${vector[$j]}"`
      done
   done
   echo ${re[@]}
}
function vectoradd(){
   vectorone=($1);
   vectortwo=($2);
   num=${#vectorone[@]};
   num=`expr $num - 1`;
   declare -a temp;
   for i in `seq 0 $num`
   do
      temp[$i]=`bc -l <<< "${vectorone[$i]}+${vectortwo[$i]}"`;
   done
   echo ${temp[@]};
}
function vectordiv(){
   vector=($1);
   scalar=$2;
   num=${#vector[@]};
   num=`expr $num - 1`;
   declare -a temp;
   for i in `seq 0 $num`
   do
      temp[$i]=`bc -l <<< "${vector[$i]}/$scalar"`;
   done
   echo ${temp[@]}
}
atom_num=5;
declare -a mass=()
#add Oxygen mass in the array#
for i in `seq 1 $( expr $atom_num \/ 5 \* 3 )`;
do
	mass+=("15.999")
done
# add Zr mass in the mass array
for i in `seq 1 $( expr $atom_num \/ 5 \* 1 )`
do
	mass+=("91.224")
done
# add oxygen mass in the mass array
for i in `seq 1 $( expr $atom_num \/ 5 \* 1 )`
do
	mass+=("137.33")
done
file_name_born="bzo.dyn1";
file_name_mod="dynmat.mold";
flag="freq (";
max_count=`grep "$flag" $file_name_born | wc -l`;
declare -A atom_i_alpha_beta;
for mu in `seq 1 $max_count`
do
   modecharge="0 0 0";
   allvector=`grepnorm $mu dynmat.mold ${mass[@]}`;
   allvector=($allvector);
   for j in `seq 1 $atom_num`
   do
      atom_i_alpha_beta=`effectcharge $j $atom_num $file_name_born`;
      vector_j="${allvector[`expr 3 \* \( $j - 1 \)`]} ${allvector[`expr 3 \* \( $j - 1 \) + 1`]} ${allvector[`expr 3 \* \( $j - 1 \) + 2`]}";
      contribution_j=`matrixtime "${atom_i_alpha_beta[@]}" "${vector_j[@]}"`;
      tick=`expr $j - 1`
      mass_sqrt=`bc -l <<< "sqrt(${mass[$tick]})"`
      contribution_j=`vectordiv "${contribution_j[@]}" $mass_sqrt`
      modecharge=`vectoradd "$modecharge" "${contribution_j[@]}"`
   done
   echo ${modecharge[@]}
done
