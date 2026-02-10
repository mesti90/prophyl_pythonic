r_container="/node8_R10/vasarhelyib/containers/stitam-prophyl-0.13.img"
treeshrink_container="/node8_R10/vasarhelyib/containers/mesti90-treeshrink.1.3.9.sif"
input_dir="tree_input"

RUN_DATING=false
RUN_DRAW=false

for arg in "$@"; do
	case "$arg" in
		--dating)
			RUN_DATING=true
			;;
		--draw)
			RUN_DRAW=true
			;;
		*)
			echo "Unknown option: $arg" >&2
			exit 1
			;;
	esac
done

function singularity_run {
        singularity run -B /node8_R10,/node8_data,/node10_R10,/home,/scratch "$@"
}

export -f singularity_run

function process_st {
	#st=ST97
	local st=$1
	local wd=work/${st}
	local gubbins_dir=${wd}/gubbins
	
	
	local gubbinstree=${gubbins_dir}/${st}.gubbins.final_tree.tre
	local gentypes=${input_dir}/${st}.tsv
	local snps=${gubbins_dir}/${st}.gubbins.filtered_polymorphic_sites.fasta
	local cpus=5
	local prophyl_dir="prophyl_custom_20260115"
	local pruned_tree=${wd}/${st}.pruned.tre
	local dropped_tips=${wd}/${st}.pruning.dropped_tips.tsv
	
	

	local treeshrink_prefix=${st}.treeshrink
	local shrinked_tree=${wd}/${treeshrink_prefix}.tre
	local rooted_prefix=${wd}/${st}.rooted
	local rooted_trees=${rooted_prefix}.trees.rds

	local dated_trees_rds=${wd}/${st}.dated_trees.rds
	local dated_trees_dir=${wd}/${st}.dated_trees
	local final_dated_tree_rds=${wd}/${st}.final_dated_tree.rds
	local final_dated_tree_nwk=work/${st}.final_dated_tree.nwk
	
	mkdir -p ${dated_trees_dir}
	mkdir -p ${wd}
	mkdir -p ${gubbins_dir}

	
	if [[ ! -s "${pruned_tree}" ]]; then
		singularity_run ${r_container} Rscript ${prophyl_dir}/bin/prune_root.R \
			--project_dir ${prophyl_dir} \
			--tree ${gubbinstree} \
			--gentypes ${gentypes} \
			--step_threshold 0.01 \
			--overall_threshold 0.01 \
			--threads ${cpus} \
			--outtree ${pruned_tree} \
			--dropped_tips ${dropped_tips}
	else
		echo "[${st}] Skipping prune_root.R: ${pruned_tree} already exists"
	fi

	if [[ ! -s "${shrinked_tree}" ]]; then
		singularity_run ${treeshrink_container}  run_treeshrink.py --tree ${pruned_tree} \
			--centroid \
			--quantiles 0.1 \
			--outprefix ${treeshrink_prefix} \
			--outdir ${wd}
	else
		echo "[${st}] Skipping run_treeshrink.py: ${shrinked_tree} already exists"
	fi

	singularity_run ${r_container} Rscript ${prophyl_dir}/bin/validate_pruning.R \
		--tree ${gubbinstree} \
		--pruned_tree ${shrinked_tree}

	singularity_run ${r_container} Rscript ${prophyl_dir}/bin/root_tree.R \
		--project_dir ${prophyl_dir} \
		--tree ${shrinked_tree} \
		--assemblies ${gentypes} \
		--root_method "all" \
		--root_topn 1 \
		--threads ${cpus} \
		--outprefix ${rooted_prefix}
		

	singularity_run ${r_container} Rscript ${prophyl_dir}/bin/date_tree.R \
		--project_dir ${prophyl_dir} \
		--trees ${rooted_trees} \
		--snps ${snps} \
		--assemblies ${gentypes} \
		--threads ${cpus} \
		--branch_dimension "snp_per_genome" \
		--clock "strict" \
		--dated_trees ${dated_trees_rds} \
		--dated_trees_dir ${dated_trees_dir} \
		--reroot false
		
	singularity_run ${r_container}  Rscript ${prophyl_dir}/bin/choose_dated_tree.R --trees ${dated_trees_rds} --out_tree_rds ${final_dated_tree_rds} --out_tree_nwk ${final_dated_tree_nwk}
}


st_list="ST105 ST22 ST30 ST398 ST45 ST59 ST5 ST8 ST97 ST9"

set -x
if $RUN_DATING; then
	export -f process_st
	for ST in ${st_list}; do
		process_st "${ST}"
	done
fi




if $RUN_DRAW; then
	SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
	
	for ST in ${st_list}; do
		singularity_run ${r_container} Rscript "${SCRIPT_DIR}/tree_drawing.R" \
			--tree "work/${ST}.final_dated_tree.nwk" \
			--meta "tree_input/${ST}.tsv" \
			--columns spatyper,Capsule.type,country,sccmec.type,PH4,PH12,PH18 \
			--out "work/${ST}.tree.pdf"
	done
fi

set +x
