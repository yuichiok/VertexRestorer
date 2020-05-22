#!/bin/bash

inputxml=VRsteer.xml

# small or large
model=l5
#model=s5

gearxml=gear_ILD_${model}_v02.xml

# yyxyev
#processID=108670
# yyxylv
#processID=108675

# eLpR
#processID=eL.pR
#processID=I108670

# type isolated lepton
isolep=electron
#isolep=muon
#isolep=tau
#isolep=muon_tau

chiral=eLpR

##===== New sample (large, new vertexing) =====##

dirin=/hsm/ilc/users/yokugawa/
processID1=default1
processID2=default2

if [ "$model" == "l5" ]; then
	if [ "$isolep" == "electron" ]; then
		dirin=/group/ilc/users/yokugawa/TTbar/l5/eLpR/semiLep/electron/slcio/
		processID1=I108668
		processID2=I108670
	else
		dirin=/group/ilc/users/yokugawa/TTbar/l5/eLpR/semiLep/muon/slcio/
		processID1=yyvlyx.eL.pR
		processID2=yyxylv.eL.pR
	fi
elif [ "$model" == "s5" ]; then
	if [ "$isolep" == "electron" ]; then
		dirin=/group/ilc/users/yokugawa/TTbar/s5/eLpR/semiLep/electron/slcio/
		processID1=I108668
		processID2=I108670
	else
		dirin=/group/ilc/users/yokugawa/TTbar/l5/eLpR/semiLep/electron/slcio/
		processID1=yyvlyx.eL.pR
		processID2=yyxylv.eL.pR
	fi
fi


# electron sample
#dirin=/hsm/ilc/users/yonamine/physics/lcfiplus/vertexing_6f_ttbar/

# muon and tau sample
#dirin=/home/ilc/yonamine/work/Yuichi/test20190328/run/slcio/

nfile=10
#nfile=12
#nfile=1
#nperjob=50

nperjob=3
./submit2 $dirin $processID1 $processID2 1 $nfile $nperjob $inputxml $gearxml $isolep $model $chiral
