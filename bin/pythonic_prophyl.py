#!/usr/bin/env python3

#TODO: known issue: gubbins puts it output in the current directory

"""
treebuilder_pipeline.py

Genome tree-building pipeline.

Usage:
    python treebuilder_pipeline.py \
        --assemblies assemblies.tsv \
        --projectID S123 \
        --threads 8 \
        [--force]
"""

import os
import csv
import re
import sys
import argparse
import shutil
import random
import pandas as pd
import numpy as np
from pathlib import Path
from spython.main import Client
from Bio import SeqIO
import pdb
import gzip
import tempfile
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from collections import defaultdict
import logging
import time
from contextlib import contextmanager
import contextlib
from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional
import subprocess

def init_logging():
	logging.basicConfig(format='%(asctime)s. %(levelname)s:%(message)s', filename=os.path.join(f"{Path(__file__).stem}_{time.strftime('%Y%m%d%H%M%S')}.log"), level=logging.INFO) #Ez ele nem kerulhet logging (msg/sysexec) parancs!!!
	logging.Formatter(fmt='%(asctime)s',datefmt="%Y")


def msg(cmd, tipus="info"):
	'''prints a message to stdout and to the log'''
	if tipus == "info":
		logging.info(cmd)
	elif tipus == "error":
		logging.error(cmd)
	timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
	print(f'{timestamp} [{tipus.upper()}] {cmd}')

def error(uzenet):
	msg(uzenet, tipus="error")
	

@dataclass
class STInfo:
	st: str
	assembly_dir: Path
	workdir: Path
	assembly_df: pd.DataFrame = field(init=False)
	
	def __post_init__(self):
		self.assembly_df = pd.read_csv(self.assemblies_file, sep="\t")
		ensure_dir(self.stdir)
		ensure_dir(self.dated_trees_dir)
	
	@cached_property
	def stdir(self) -> Path:
		return self.workdir / f"ST{self.st}"
	
	@cached_property
	def assemblies_file(self) -> Path:
		return self.assembly_dir / f"ST{self.st}.tsv"
		
	@cached_property
	def refgenome(self) -> Path:
		return self.stdir / f"ST{self.st}.reference.fna"
	
	@cached_property
	def snippy_dir(self) -> Path:
		return self.stdir / "snippy"
		
	@cached_property
	def gubbins_dir(self) -> Path:
		return self.stdir / "gubbins"
		
	@cached_property
	def snippy_out(self) -> Path:
		return self.stdir / f"ST{self.st}.snippy_out.fna"
	
	@cached_property
	def snippy_dedup(self) -> Path:
		return self.stdir / f"ST{self.st}.snippy_out.nodup.fna"
	
	@cached_property
	def reference_strain_file(self) -> Path:
		return self.stdir / "reference_strain.txt"
	
	@cached_property
	def gubbins_prefix(self) -> str:
		return f"ST{self.st}.gubbins"
	
	@cached_property
	def gubbins_tree(self) -> Path:
		return self.gubbins_dir / f"{self.gubbins_prefix}.final_tree.tre"
	
	@cached_property
	def gubbins_fasta(self) -> Path:
		return self.gubbins_dir / f"{self.gubbins_prefix}.filtered_polymorphic_sites.fasta"
	
	@cached_property
	def treeshrink(self) -> Path:
		return self.stdir / f"{self.treeshrink_prefix}.tre"
	
	@cached_property
	def treeshrink_prefix(self):
		return f"ST{self.st}.treeshrink"
	
	@cached_property
	def rooted_shrink(self):
		return self.stdir / f"{self.treeshrink_prefix}.treeshrink.rooted.nwk"
	
	@cached_property
	def tree_prune(self):
		return self.stdir / f"ST{self.st}.treepruner.nwk"
		
	@cached_property
	def pruned_tree(self) -> Path:
		return self.stdir / f"ST{self.st}.pruned.tre"

	@cached_property
	def dropped_tips(self) -> Path:
		return self.stdir / f"ST{self.st}.pruning.dropped_tips.tsv"

	@cached_property
	def rooted_prefix(self) -> Path:
		return self.stdir / f"ST{self.st}.rooted"

	@cached_property
	def rooted_trees(self) -> Path:
		return Path(f"{self.rooted_prefix}.trees.rds")

	@cached_property
	def dated_trees_rds(self) -> Path:
		return self.stdir / f"ST{self.st}.dated_trees.rds"

	@cached_property
	def dated_trees_dir(self) -> Path:
		return self.stdir / f"ST{self.st}.dated_trees"

	@cached_property
	def final_dated_tree_rds(self) -> Path:
		return self.stdir / f"ST{self.st}.final_dated_tree.rds"

	@cached_property
	def final_dated_tree_nwk(self) -> Path:
		return self.workdir / f"ST{self.st}.final_dated_tree.nwk"
	
	@cached_property
	def tree_figure(self) -> Path:
		return self.workdir / f'ST{self.st}.tree.pdf'

# ============================================================
# Utilities
# ============================================================


container_names = {
	"snippy": "/home/vasarhelyib/containers/staphb-snippy-4.6.0-SC2.img",
	"seqkit": "/home/vasarhelyib/containers/staphb-seqkit.2.10.0.sif",
	"gubbins": "/home/vasarhelyib/containers/staphb-gubbins.3.3.5.sif",
	"r_container": "/home/vasarhelyib/containers/stitam-prophyl-0.13.img",
	"treeshrink": "/home/vasarhelyib/containers/mesti90-treeshrink.1.3.9.sif"
}

containers = {name: Client.instance(img) for name, img in container_names.items()}


def ensure_dir(path):
	Path(path).mkdir(parents=True, exist_ok=True)

def skip_if_done(path, force=False, quiet=False):
	path = Path(path)
	if path.exists() and not force:
		if not quiet:
			msg(f"[SKIP] {path} already exists.")
		return True
	return False


def run_container(inst_name, command, quiet=False):
	tmpdir = Path(tempfile.mkdtemp(prefix="spython_tmp_"))
	scratch_dir = tmpdir / "scratch" / "tmp"
	scratch_dir.mkdir(parents=True, exist_ok=True)
	bind_options = ["/node8_R10", "/node10_R10", str(tmpdir), f"{scratch_dir}:/scratch/tmp"]

	if not quiet:
		msg(f"[RUN] {container_names[inst_name]} :: {command}")

	try:
		result = Client.execute(containers[inst_name], command, bind=bind_options)
	except Exception as e:
		# Optionally log the error instead of silently passing
		print(f"[ERROR] Container execution failed: {e}")
		result = None
	finally:
		shutil.rmtree(tmpdir, ignore_errors=True)

	return result


def get_args():
	"""
	Parse command-line arguments and return them as a Namespace.
	If `args` is provided, it should be a list of strings (for testing purposes).
	"""
	parser = argparse.ArgumentParser(description="Genome treebuilder pipeline")
	parser.add_argument("--st_file", default="config/st_list.20260105.txt",help="Text file with one ST per line")
	parser.add_argument("--typing_table", default="data/sample_types.inlab_samples.phage_sensitivity.20260128.csv", help="Table containing all assembly typings, and phage sensitivities")
	parser.add_argument("--assembly_dir", default="tree_input", help="Directory storing assembly tables")
	parser.add_argument("--all_assemblies", default="data/all_assemblies.tsv", help="Table containing all assembly paths")
	parser.add_argument("-wd","--workdir",default="work")
	parser.add_argument("-n","--cpu", type=int, default=20, help="Number of threads")
	parser.add_argument("--subthreads", type=int, default=4, help="Number of subthreads in a thread")
	parser.add_argument("--force", action="store_true", help="Overwrite existing files")
	
	parser.add_argument("--snippy", action="store_true", help="Run snippy")
	parser.add_argument("--gubbins", action="store_true", help="Run gubbins")
	parser.add_argument("--date_root", action="store_true", help="Run dating and rooting")
	parser.add_argument("--scriptdir", default="/node10_R10/vasarhelyib/Staphylococcus_aureus/prophyl_pythonic/bin")
	parser.add_argument("--draw", action="store_true", help="Draw trees")
	parser.add_argument("--heatmap_vars", default="spatyper,Capsule.type,country,sccmec.type,PH4,PH12,PH18")
	parser.add_argument("--all", action="store_true", help="Run all parts")
	
	args = parser.parse_args()
	if args.all:
		args.snippy = True
		args.gubbins = True
		args.date_root = True
	args.assembly_dir = Path(args.assembly_dir)
	args.workdir = Path(args.workdir)
	args.scriptdir = Path(args.scriptdir)
	args.localdir = ".local"
	return args

# ============================================================
# Step 1: Read Assemblies
# ============================================================

def read_st_file(st_file):
	"""
	Read a list of STs from a file.
	Ignores empty lines and lines starting with '#'.
	"""
	try:
		msg(f"Reading ST file: {st_file}")
		if not os.path.isfile(st_file):
			raise FileNotFoundError(f"{st_file} is missing")
		with open(st_file) as f:
			sts = {line.strip() for line in f if line.strip() and not line.startswith("#")}
		if not sts:
			raise ValueError("ST file is empty after filtering")
		msg(f"Loaded {len(sts)} STs")
		return sts

	except Exception as e:
		msg(f"Failed to read ST file: {e}", "error")
		raise

def read_assembly_paths(assembly_table):
	"""
	Read assembly paths table and return {strain: assembly_path}.
	"""
	try:
		msg(f"Reading assembly table: {assembly_table}")
		df = pd.read_csv(assembly_table, sep="\t", usecols=["strain", "assembly_path"], dtype=str)
		if df.empty:
			raise ValueError("Assembly table is empty")
		assembly_dict = dict(zip(df["strain"], df["assembly_path"]))
		msg(f"Loaded {len(assembly_dict)} assembly paths")
		return assembly_dict
	except Exception as e:
		error(f"Failed to read assembly paths: {e}")
		raise

def read_typing_table(typing_table_file):
	"""
	Read typing table into a pandas DataFrame.
	"""
	try:
		msg(f"Reading typing table: {typing_table_file}")
		df = pd.read_csv(typing_table_file, sep="\t", dtype=str).fillna("")
		if df.empty:
			raise ValueError("Typing table is empty")
		msg(f"Typing table rows: {len(df)}")
		return df
	except Exception as e:
		error(f"Failed to read typing table: {e}")
		raise

def write_assembly_files(typing_table_df, assembly_dict, st_list, assembly_dir):
	"""
	Write one TSV file per ST containing strain and assembly_path.
	Only existing assembly files are written; missing paths are logged as errors.
	"""
	try:
		msg("Writing assembly files per ST")
		
		#Collect data from typing table - common for all STs
		filtered_df = typing_table_df[typing_table_df["mlst"].isin(st_list)]
		typing_dict = filtered_df.set_index('Genome').to_dict('index')
		typing_cols = list(filtered_df.columns.drop('Genome'))
		for col in ['genome_size','longest_contig']:
			if col not in typing_cols:
				typing_cols.append("longest_contig")
	
		#Create groups by ST
		grouped = filtered_df.groupby("mlst")
		
		#Track duplicate assembly paths
		duplicates = [["ST","representative","duplicate strains"]]
		
		for st in st_list:
			df = grouped.get_group(st) if st in grouped.groups else pd.DataFrame(columns=typing_table_df.columns) # create a fall back - if the given ST is missing (it's not possible, as st_list contains the most prevalent STs)
			
			#We want to remove (accidentally) duplicate identifiers
			seen_strains = set()
			#For all strains gather the related assembly paths
			assemblies_for_st = defaultdict(list)
			
			for strain in df["Genome"]:
				if strain in seen_strains:
					continue
				#Add the assembly path for the strain
				assembly_path = assembly_dict.get(strain)
				assembly_path = Path(assembly_path) if assembly_path else None
				if assembly_path and assembly_path.is_file():
					assemblies_for_st[assembly_path].append(strain)
					seen_strains.add(strain)
				else:
					msg(f"ERROR: Missing assembly for strain '{strain}': {assembly_path}","error")
			
			outrows = []
			
			#For all assembly path select the first strain, and add also the typing details
			for path, strains in assemblies_for_st.items():
				keep = next((s for s in strains if s.startswith("SP")), strains[0])  # pick SP* if exists, else first
				if len(strains) > 1:
					duplicates.append([st, keep, ",".join(s for s in strains if s != keep)])
				
				td = typing_dict.get(keep, {})
				genome_size = td.get("genome_size")
				longest_contig = td.get("longest_contig")
				if genome_size is None or longest_contig is None:
					try:
						genome_size, longest_contig = get_genome_metrics(path)
					except Exception as e:
						msg(f"ERROR computing metrics for {keep}: {e}", "error")
						genome_size, longest_contig = None, None

				row = {
					"strain": keep,
					"assembly_path": path,
					**{col: td.get(col) for col in typing_cols if col not in ["genome_size", "longest_contig"]},
					"genome_size": genome_size,
					"longest_contig": longest_contig
				}
				outrows.append(row)
			
			outfile = assembly_dir / f"ST{st}.tsv"
			df = (
				pd.DataFrame(outrows, columns=["strain", "assembly_path"] + typing_cols)
				.assign(
						_key1=lambda x: np.where(x["strain"].str.startswith("SPA"), 0, 1),
						_key2=lambda x: x["strain"].str.extract(r"^SPA(\d+)").astype(float).fillna(0),
						_key3=lambda x: x["strain"]
					)
				.sort_values(by=["_key1","_key2","_key3"])
				.drop(columns=["_key1","_key2","_key3"])
				.reset_index(drop=True)
			)
			df.to_csv(outfile, sep="\t", index=False)
			msg(f"{outfile} written ({len(assemblies_for_st)} assemblies)")
		
		with (assembly_dir / "duplicates.tsv").open("w") as g:
			wtr = csv.writer(g, delimiter="\t")
			wtr.writerows(duplicates)
		
		msg("All ST assembly files written")

	except Exception as e:
		error(f"Failed while writing assembly files: {e}")
		raise

def create_assembly_tables(st_list, args):
	"""
	Main workflow to generate per-ST assembly tables.
	"""
	marker_file = Path(args.localdir) / "assembly_tables_ready"
	if marker_file.exists() and not args.force:
		msg(f"Assembly tables already ready ({marker_file}). Skipping generation.")
		return
	try:
		msg("Starting assembly table creation")
		assembly_dict = read_assembly_paths(args.all_assemblies)
		typing_table_df = read_typing_table(args.typing_table)
		ensure_dir(args.assembly_dir)
		write_assembly_files(typing_table_df, assembly_dict, st_list, args.assembly_dir)
		msg("Pipeline finished successfully")
		timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
		marker_file.write_text(f'{timestamp} Assembly tables generated successfully.\n')

	except Exception as e:
		logging.critical(f"Pipeline failed: {e}")
		raise





# ============================================================
# Step 3: Choose Reference Genome
# ============================================================

def get_genome_metrics(fasta_path):
	"""Return genome size and longest contig length for a FASTA (plain or gzipped)."""
	if not os.path.exists(fasta_path):
		raise FileNotFoundError(f"Missing FASTA: {fasta_path}")

	# Handle gzipped input
	open_func = gzip.open if fasta_path.endswith(".gz") else open
	with open_func(fasta_path, "rt") as handle:
		sizes = [len(rec.seq) for rec in SeqIO.parse(handle, "fasta")]

	if not sizes:
		raise ValueError(f"No contigs found in {fasta_path}")
	return sum(sizes), max(sizes)


def choose_reference_genome(info, args):
	if skip_if_done(info.refgenome, args.force):
		return

	msg(f"ST{info.st}: Selecting reference genome...")
	# Ensure copy to avoid pandas view warnings
	df = info.assembly_df.copy()

	# If there are GCF/GCA genomes, keep only those
	mask = df["strain"].str.contains(r"GC[AF]_", regex=True)
	if mask.any():
		df = df[mask]
		msg(f"ST{info.st}: Kept only GCF/GCA assemblies ({len(df)} remain)")


	# Drop genome size outliers (keep middle 95%)
	lower, upper = np.percentile(df["genome_size"], [2.5, 97.5])
	df = df[(df["genome_size"] >= lower) & (df["genome_size"] <= upper)]
	msg(f"ST{info.st}: Kept {len(df)} assemblies in genome size 95% range.")

	# Select by longest contig, earliest date, random tie-break
	df = df.sort_values(by=["longest_contig", "Date"], ascending=[False, True])
	top = df.iloc[0]

	# Handle ties by contig length and date
	same_len = df[df["longest_contig"] == top["longest_contig"]]
	if len(same_len) > 1:
		top = same_len.sort_values("collection_date").iloc[0]
		tie_same_date = same_len[same_len["collection_date"] == top["collection_date"]]
		if len(tie_same_date) > 1:
			top = same_len.sample(1).iloc[0]

	chosen_path = Path(top["assembly_path"])

	# gunzip if necessary when copying reference genome
	msg(f"ST{info.st}: Copying reference genome to {info.refgenome}")
	open_func = gzip.open if chosen_path.suffix == ".gz" else open
	
	with open_func(chosen_path, "rb") as f_in, open(info.refgenome, "wb") as f_out:
		shutil.copyfileobj(f_in, f_out)

	msg(f"ST{info.st}: Selected reference genome: {top['strain']} → {info.refgenome}")
	with open(info.stdir / "reference_strain.txt", "w") as g:
		g.write(top['strain'])



# ============================================================
# Step 4: Run Snippy per assembly
# ============================================================

def run_snippy_single(assembly_id, assembly_path, snippy_dir, ref_genome):
	snippy_dir = Path(snippy_dir)
	ref_genome = Path(ref_genome)
	assembly_path = Path(assembly_path)
	
	outfile = snippy_dir / f"{assembly_id}.snippy.fna"
	tmp_ctg = snippy_dir / f"{assembly_id}.tmp.fna"
	
	if skip_if_done(outfile, force=False, quiet=True):
		return

	# Decompress or copy
	try:
		if assembly_path.endswith(".gz"):
			with gzip.open(assembly_path, "rb") as f_in, open(tmp_ctg, "wb") as f_out:
				shutil.copyfileobj(f_in, f_out)
		else:
			shutil.copy(assembly_path, tmp_ctg)

		cmd = f"snippy --outdir {snippy_dir/assembly_id} --ref {ref_genome} --ctgs {tmp_ctg} --force"
		run_container("snippy", cmd, quiet=True)

		# look for possible outputs
		snps_path = snippy_dir / assembly_id / "snps.consensus.subs.fa"
		if os.path.exists(snps_path):
			shutil.move(snps_path, outfile)

	except Exception as e:
		print(f"[ERROR] Snippy failed for {assembly_id}: {e}")

	finally:
		tmpdir = snippy_dir / assembly_id
		#if os.path.exists(tmpdir):
		#	shutil.rmtree(tmpdir, ignore_errors=True)
		#if os.path.exists(tmp_ctg):
		#	os.remove(tmp_ctg)


def run_snippy_all(df, snippy_dir, ref_genome, threads, st):
	msg(f"[ST{st}] Running Snippy on {len(df)} assemblies (threads={threads})...")
	if threads == 1:
		# sequential
		for _, row in df.iterrows():
			run_snippy_single(row["strain"], row["assembly_path"], snippy_dir, ref_genome)
	else:
		# multithreaded
		with ThreadPoolExecutor(max_workers=threads) as executor:
			futures = {executor.submit(run_snippy_single, row["strain"], row["assembly_path"], snippy_dir, ref_genome): row["strain"] for _, row in df.iterrows()}
			for future in as_completed(futures):
				future.result() 
	# final summary
	msg(f"[ST{st}] Snippy finished")


def concat_snippy(info): 
	"""
	Concatenate the longest contig from each per-assembly Snippy file into a single FASTA.
	The output FASTA format is simplified: one line per sequence.
	"""
	msg(f"[ST{info.st}] concatenating snippy results")
	snippy_files = []
	for _, row in info.assembly_df.iterrows():
		assembly_id = row["strain"]
		snippy_file = info.snippy_dir / f"{assembly_id}.snippy.fna"
		if snippy_file.exists():
			snippy_files.append((assembly_id, snippy_file))
		else:
			msg(f"[WARN] [ST{info.st}] Snippy output not found for {assembly_id}, skipping")

	if not snippy_files:
		msg(f"[WARN] [ST{info.st}] No Snippy files found to concatenate.")
		return
	
	msg(f"[ST{info.st}] {len(snippy_files)} are to be concatenated")

	try:
		with open(info.snippy_out, "w") as out_f:
			for assembly_id, f in snippy_files:
				# Find the longest contig
				contigs = list(SeqIO.parse(f, "fasta"))
				if not contigs:
					msg(f"[WARN] No contigs in {f}, skipping")
					continue
				longest_contig = max(contigs, key=lambda x: len(x.seq))
				# Write in single-line format
				seq_line = str(longest_contig.seq)
				out_f.write(f">{assembly_id}\n{seq_line}\n")
	except Exception as e:
		msg(f"Failed to write concatenated file: {e}", "error")
		return
	msg(f"Concatenated longest contigs → {info.snippy_out} ({len(snippy_files)} assemblies)")
	# Only remove individual files after successful concatenation
	for _, f in snippy_files:
		try:
			f.unlink()
		except Exception as e:
			msg(f"[WARN] Could not remove {f}: {e}")







# ============================================================
# Step 6: Remove duplicates using seqkit
# ============================================================

def remove_duplicates(info, args):
	if skip_if_done(info.snippy_dedup, args.force):
		return info.snippy_dedup
	cmd = f"seqkit rmdup -s {info.snippy_out} -D {info.snippy_out.parent}/duplicates.txt -o {info.snippy_dedup} -j {args.subthreads}"
	run_container("seqkit", cmd)
	print(f"[OK] Duplicates removed → {info.snippy_dedup}")


# ============================================================
# Step 7: Run Gubbins
# ============================================================

def run_gubbins(info, args): 
	ensure_dir(info.gubbins_dir)
	required_files = [info.gubbins_tree, info.gubbins_fasta]
	
	def required_exist(directory):
		return all((directory / fpath.name ).exists() for fpath in required_files)
	
	#If the final results exist, then skip processing
	if required_exist(info.gubbins_dir):
		msg(f"[SKIP] [ST{info.st}] gubbins is already done")
		return
	
	cwd = Path.cwd()
	if not required_exist(cwd):
		cmd = f"run_gubbins.py --model-fitter raxmlng --tree-builder fasttree --threads {args.subthreads} --prefix {info.gubbins_prefix} --iterations 10 {Path(info.snippy_dedup).resolve()}"
		run_container("gubbins", cmd)
	
	##Move gubbins output files to workdir!!!
	output_suffices = [
		"branch_base_reconstruction.embl",
		"filtered_polymorphic_sites.phylip",
		"log",
		"node_labelled.final_tree.tre",
		"per_branch_statistics.csv",
		"recombination_predictions.embl",
		"recombination_predictions.gff",
		"summary_of_snp_distribution.vcf"
	]
	
	for suffix in output_suffices:
		name = f"{info.gubbins_prefix}.{suffix}"
		src = cwd / name
		if src.exists():
			shutil.move(src, info.gubbins_dir / name)
	for fpath in required_files:
		src = cwd / fpath.name
		if src.exists():
			shutil.move(cwd / fpath.name, fpath)

	msg(f"[ST{info.st}] Gubbins finished")



# ============================================================
# Step 8–11: Tree operations (R-based)
# ============================================================




# ============================================================
# Main
# ============================================================


def run_parallel(st_info_list, func, max_workers=4, **kwargs):
	"""
	Run a stage function in parallel over a list of STInfo objects.
	func: a function that takes (STInfo, **kwargs)
	"""
	results = []
	with ThreadPoolExecutor(max_workers=max_workers) as executor:
		futures = {executor.submit(func, info, **kwargs): info.st for info in st_info_list}
		for future in as_completed(futures):
			st = futures[future]
			results.append(future.result())
	return results




def collect_all_snippy_tasks(st_info_list):
	for info in st_info_list:
		ensure_dir(info.snippy_dir)
		for _, row in info.assembly_df.iterrows():
			yield (info, row["strain"], row["assembly_path"])

def run_snippy_all_parallel(st_info_list, args):
	tasks = list(collect_all_snippy_tasks(st_info_list))

	msg(f"[INFO] Running Snippy on {len(tasks)} assemblies (threads={args.cpu})")

	with ThreadPoolExecutor(max_workers=args.cpu) as executor:
		futures = [
			executor.submit(run_snippy_single, 	strain, assembly_path, info.snippy_dir, info.refgenome,)
			for info, strain, assembly_path in tasks
		]

		for future in as_completed(futures):
			future.result()

def finalize_snippy_st(info, args):
	if not skip_if_done(info.snippy_out, args.force):
		concat_snippy(info)

	remove_duplicates(info, args)

def run_R(script, args_line, args, container="r_container"):
	print(f"[INFO] Running R step: {script}")
	subprocess.run(f"singularity run --bind /node8_R10,/node10_R10 {container_names[container]} Rscript {args.scriptdir / script} {args_line}", shell=True)



def prune_tree(info, args):
	if not info.pruned_tree.exists() or info.pruned_tree.stat().st_size == 0:
		cmd=f"--project_dir {args.scriptdir} --tree {info.gubbins_tree} --gentypes {info.assemblies_file} --step_threshold 0.01 --overall_threshold 0.01 --threads {args.subthreads} --outtree {info.pruned_tree} --dropped_tips {info.dropped_tips}"
		run_R("prune_root.R", cmd, args=args)
	else:
		print(f"[{info.st}] Skipping prune_root.R: {info.pruned_tree} exists")


def shrink_tree(info, args):
	# ---- TreeShrink ----
	if not info.treeshrink.exists() or info.treeshrink.stat().st_size == 0:
		cmd = (
			f"run_treeshrink.py "
			f"--tree {info.pruned_tree} "
			f"--centroid "
			f"--quantiles 0.1 "
			f"--outprefix {info.treeshrink_prefix} "
			f"--outdir {info.stdir}"
		)
		run_container(treeshrink_container, cmd)
	else:
		print(f"[{info.st}] Skipping TreeShrink: {info.treeshrink} exists")

def validate_pruning(info, args):
	run_R("validate_pruning.R", f"--tree {info.gubbins_tree} --pruned_tree {info.treeshrink}", args=args)


def root_tree(info, args):
	if info.rooted_trees.exists():
		msg(f"[{info.st}] Skip rooting tree, as {info.rooted_trees} exists")
		return
	run_R(
		"root_tree.R",
		f"--project_dir {args.scriptdir} "
		f"--tree {info.treeshrink} "
		f"--assemblies {info.assemblies_file} "
		f"--outprefix {info.rooted_prefix} "
		f"--root_method all "
		f"--root_topn 1 "
		f"--threads {args.subthreads}",
		args=args,
	)

def date_tree(info, args):
	if info.dated_trees_rds.exists():
		msg(f"[{info.st}] Skip dating tree, as {info.dated_trees_rds} exists")
		return
	run_R(
		"date_tree.R",
		f"--project_dir {args.scriptdir} "
		f"--trees {info.rooted_trees} "
		f"--snps {info.gubbins_fasta} "
		f"--assemblies {info.assemblies_file} "
		f"--threads {args.subthreads} "
		f"--branch_dimension snp_per_genome "
		f"--clock strict "
		f"--reroot false "
		f"--dated_trees {info.dated_trees_rds} "
		f"--dated_trees_dir {info.dated_trees_dir}",
		args=args,
	)

def choose_dated_tree(info, args):
	if info.final_dated_tree_rds.exists() and info.final_dated_tree_nwk.exists():
		msg(f"[{info.st}] Skip dating tree, as {info.final_dated_tree_rds} and {info.final_dated_tree_nwk} exist")
		return
	run_R(
		"choose_dated_tree.R",
		f"--trees {info.dated_trees_rds} "
		f"--out_tree_rds {info.final_dated_tree_rds} "
		f"--out_tree_nwk {info.final_dated_tree_nwk}",
		args=args,
	)


def draw_tree(info, args):
	if info.tree_figure.exists():
		msg(f"[{info.st}] Skip drawing tree, as {info.tree_figure} exists")
		return
	run_R(
		"tree_drawing.R",
		f"--tree {info.final_dated_tree_nwk} "
		f"--meta {info.assemblies_file} "
		f"--columns {args.heatmap_vars} "
		f"--out {info.tree_figure} "
		f"--scriptdir {args.scriptdir}",
		args=args,
	)


def date_and_root(info, args):
	# ---- prune_root.R ----
	prune_tree(info, args)
	shrink_tree(info, args)
	validate_pruning(info, args)
	root_tree(info, args)
	date_tree(info, args)
	choose_dated_tree(info, args)


def main():
	init_logging()
	args = get_args()
	# Setup & prepare assemblies
	ensure_dir(args.localdir)
	ensure_dir(args.workdir)
	
	#Get ST list
	st_list = read_st_file(args.st_file)
	create_assembly_tables(st_list, args)
	
	st_info_list = [STInfo(st = st, assembly_dir = args.assembly_dir, workdir = args.workdir) for st in st_list]
	
	#Snippy
	if args.snippy:
		msg("====Snippy====")
		run_parallel(st_info_list, choose_reference_genome, max_workers=args.cpu, args=args)
		run_snippy_all_parallel(st_info_list, args)
		run_parallel(st_info_list, finalize_snippy_st, max_workers=args.cpu, args=args)
	
	#Gubbins
	if args.gubbins:
		msg("====Gubbins====")
		run_parallel(st_info_list, run_gubbins,  args=args)
	
	#Post-process
	if args.date_root:
		msg("====Dating and rooting trees====")
		run_parallel(st_info_list, date_and_root, max_workers=args.cpu, args=args)
	
	if args.draw:
		msg("====Draw trees====")
		run_parallel(st_info_list, draw_tree, args=args)


if __name__ == "__main__":
	main()


