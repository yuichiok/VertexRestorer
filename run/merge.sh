#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	echo "Please enter isolated lepton type."
	exit 1
fi

model=$1
isolep=$2
chiral=$3
#dir=/hsm/ilc/users/yokugawa/preset_N_run/${model}/semiLep/${chiral}/${isolep}
dir=/group/ilc/users/yokugawa/TTbar/${model}/${chiral}/semiLep/${isolep}

if [ ! -d "$dir" ]; then
	echo "The directory does not exist."
	exit 1
else
	echo ${dir}
	#hadd -f ${dir}/QQbarProcessor_out/root_merge/${chiral}_${isolep}_QQbar.root ${dir}/QQbarProcessor_out/*.root
	#hadd -f ${dir}/TrashRecoProcessor_out/before_vtx_recovery/root_merge/${chiral}_${isolep}_TRP_before_D0.root ${dir}/TrashRecoProcessor_out/before_vtx_recovery/*.root
	hadd -f ${dir}/TrashRecoProcessor_out/after_vtx_recovery/root_merge/${chiral}_${isolep}_TRP_after_test.root ${dir}/TrashRecoProcessor_out/after_vtx_recovery/*.root
	#hadd -f ${dir}/particle_tagger_out/root_merge/${chiral}_${isolep}_PT_D0.root ${dir}/particle_tagger_out/*.root
	#hadd -f ${dir}/truth_vertex_finder_out/root_merge/${chiral}_${isolep}_TVF_D0.root ${dir}/truth_vertex_finder_out/*.root
	hadd -f ${dir}/vertex_restorer_out/root_merge/${chiral}_${isolep}_VR_test.root ${dir}/vertex_restorer_out/*.root
fi
