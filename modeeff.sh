#!/bin/bash
function grepnorm(){
	count=$1;
	atom_num=$2;
	file_name=$3;
    mass=$4;
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
    ####################### send result out ########################################################
    declare -a re;
    tick=`expr $atom_num - 1`;
    tickone=`expr $tick \* 3`;
    ticktwo=`expr $tickone + 1`;
    tickthree=`expr $tickone + 2`;
    re+=(${allvector[$tickone]});
    re+=(${allvector[$ticktwo]});
    re+=(${allvector[$tickthree]});
    rm temp1.txt
    echo ${re[@]}
}
function effecharge(){
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
atom_num=5;
declare -a mass=();
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
declare -A modecharge;
for mu in `seq 1 $max_count`
do
    for alpha in `seq 0 2`
    do
 	    modecharge[$alpha]="0.0";
 	    #please see the formula in Bennet's paper
 	    for i in `seq 1 $atom_num`
 	    do
 		    z_i_alpha_beta=`effecharge $i $atom_num $file_name_born | sed -n "$( expr 1 + $alpha )p"`;
            z_i_alpha_beta=($z_i_alpha_beta);
 		    for beta in `seq 0 2`
 		    do
 			    z_i_alpha_beta=${z_i_alpha_beta[$beta]};
 		    	a_mu_i_beta=`grepnorm $mu $atom_num $file_name_mod $mass`;
 		    	a_mu_i_beta=($a_mu_i_beta);
 		    	a_mu_i_beta=${a_mu_i_beta[$beta]};
 		    	add=`bc -l <<< "$z_i_alpha_beta*$a_mu_i_beta/(sqrt(${mass[$( expr $i - 1 )]}))"`
 		    	modecharge[$alpha]=`bc -l <<< "${modecharge[$alpha]}+$add" | awk '{printf "%2.14f", $0}'`;
 		    done
 	    done
    done
    echo ${modecharge[@]}
done
