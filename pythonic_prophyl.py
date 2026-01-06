#!/usr/bin/env python3
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
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

# ============================================================
# Utilities
# ============================================================


container_names = {
	"snippy": "/home/vasarhelyib/containers/staphb-snippy-4.6.0-SC2.img",
	"seqkit": "/home/vasarhelyib/containers/staphb-seqkit.2.10.0.sif",
	"gubbins": "/home/vasarhelyib/containers/staphb-gubbins.3.3.5.sif",
	"r_container": "/home/vasarhelyib/containers/stitam-prophyl-0.13.img"
}

containers = {name: Client.instance(img) for name, img in container_names.items()}


def ensure_dir(path: Path):
	path.mkdir(parents=True, exist_ok=True)

def skip_if_done(path: Path, force: bool, quiet=False) -> bool:
	if path.exists() and not force:
		if not quiet:
			(f"[SKIP] {path} already exists.")
		return True
	return False


def run_container(inst_name, command: str, quiet=False):
	tmpdir = Path(tempfile.mkdtemp(prefix="spython_tmp_"))
	scratch_dir = tmpdir / "scratch" / "tmp"
	scratch_dir.mkdir(parents=True, exist_ok=True)
	bind_options = ["/node8_R10", "/node10_R10", str(tmpdir), f"{scratch_dir}:/scratch/tmp"]

	if not quiet:
		print(f"[RUN] {container_names[inst_name]} :: {command}")

	try:
		result = Client.execute(containers[inst_name], command, bind=bind_options)
	except Exception as e:
		# Optionally log the error instead of silently passing
		print(f"[ERROR] Container execution failed: {e}")
		result = None
	finally:
		shutil.rmtree(tmpdir, ignore_errors=True)

	return result


def get_args(args=None):
	"""
	Parse command-line arguments and return them as a Namespace.
	If `args` is provided, it should be a list of strings (for testing purposes).
	"""
	parser = argparse.ArgumentParser(description="Genome treebuilder pipeline")
	parser.add_argument("--assemblies", default="input/ST117.tsv", help="TSV file listing assemblies")
	parser.add_argument("--projectID", default="ST117", help="Project ID")
	parser.add_argument("--threads", type=int, default=30, help="Number of threads")
	parser.add_argument("--force", action="store_true", help="Overwrite existing files")
	parser.add_argument("--root_method", default="midpoint", help="Tree rooting method")
	parser.add_argument("--root_topn", type=int, default=5, help="Top N candidates for rooting")
	parser.add_argument("--clock", default="strict", help="Molecular clock model")

	return parser.parse_args(args)

# ============================================================
# Step 1: Read Assemblies
# ============================================================

def read_assemblies(assemblies_file: Path):
	df = pd.read_csv(assemblies_file, sep="\t", usecols=["assembly", "assembly_path", "collection_date"])
	print(f"[INFO] Loaded {len(df)} assemblies.")
	return df


# ============================================================
# Step 2: Validate Assemblies
# ============================================================

def validate_assemblies(df, workdir: Path, force=False):
	outfile = workdir / "validated.info"
	if os.path.isfile(outfile):
		with open(outfile) as f:
			return f.read().strip() == "1"
	success = True
	
	print("[INFO] Validating assemblies...")
	dupes = df[df["assembly"].duplicated(keep=False)]
	if not dupes.empty:
		print("[ERROR] Duplicate assembly names found:\n")
		print(dupes["assembly"].to_string(index=False))
		success = False

	for p in df["assembly_path"]:
		if not os.path.exists(p):
			print(f"[ERROR] Missing assembly file: {p}")
			success = False
	if success:
		print(f"[OK] Validation complete, {outfile} is saved")
		with open(outfile, "w") as g:
			g.write("1")
	
	return success


# ============================================================
# Step 3: Choose Reference Genome
# ============================================================

def get_genome_metrics(fasta_path: Path):
	"""Return genome size and longest contig length for a FASTA (plain or gzipped)."""
	if not fasta_path.exists():
		raise FileNotFoundError(f"Missing FASTA: {fasta_path}")

	# Handle gzipped input
	if str(fasta_path).endswith(".gz"):
		handle = gzip.open(fasta_path, "rt")
	else:
		handle = open(fasta_path, "r")

	with handle:
		sizes = [len(rec.seq) for rec in SeqIO.parse(handle, "fasta")]

	if not sizes:
		raise ValueError(f"No contigs found in {fasta_path}")

	genome_size = sum(sizes)
	longest_contig = max(sizes)
	return genome_size, longest_contig




def choose_reference_genome(df, workdir: Path, projectID: str, force=False):
	outfile = workdir / f"{projectID}.reference_genome.fna"
	if skip_if_done(outfile, force):
		print(f"[INFO] Reference genome already exists → {outfile}")
		return outfile

	print("[INFO] Selecting reference genome...")

	# Ensure copy to avoid pandas view warnings
	df = df.copy()

	# 1. If there are GCF/GCA genomes, keep only those
	if df["assembly"].str.contains(r"GC[AF]_", regex=True).any():
		df = df[df["assembly"].str.contains(r"GC[AF]_", regex=True)]
		print(f"[INFO] Kept only GCF/GCA assemblies ({len(df)} remain)")

	# 2. Compute genome size and longest contig (batch assignment)
	genome_sizes = []
	longest_contigs = []

	for path in df["assembly_path"]:
		size, contig = get_genome_metrics(Path(path))
		genome_sizes.append(size)
		longest_contigs.append(contig)

	df["genome_size"] = genome_sizes
	df["longest_contig"] = longest_contigs

	# 3. Drop genome size outliers (keep middle 95%)
	lower, upper = np.percentile(df["genome_size"], [2.5, 97.5])
	df = df[(df["genome_size"] >= lower) & (df["genome_size"] <= upper)]
	print(f"[INFO] Kept {len(df)} assemblies in genome size 95% range.")

	# 4. Select by longest contig, earliest date, random tie-break
	df = df.sort_values(by=["longest_contig", "collection_date"], ascending=[False, True])
	top = df.iloc[0]

	# Handle ties by contig length and date
	if len(df[df["longest_contig"] == top["longest_contig"]]) > 1:
		same_len = df[df["longest_contig"] == top["longest_contig"]]
		top = same_len.sort_values("collection_date").iloc[0]
		if len(same_len[same_len["collection_date"] == top["collection_date"]]) > 1:
			top = same_len.sample(1).iloc[0]

	chosen_path = Path(top["assembly_path"])

	# gunzip if necessary when copying reference genome
	print(f"[INFO] Copying reference genome → {outfile}")
	if str(chosen_path).endswith(".gz"):
		with gzip.open(chosen_path, "rb") as f_in, open(outfile, "wb") as f_out:
			shutil.copyfileobj(f_in, f_out)
	else:
		shutil.copy(chosen_path, outfile)

	print(f"[OK] Selected reference genome: {top['assembly']} → {outfile}")
	return outfile



# ============================================================
# Step 4: Run Snippy per assembly
# ============================================================

def run_snippy_single(assembly_id: str, assembly_path: Path, snippy_dir: Path, ref_genome: Path, force=False):
	outfile = snippy_dir / f"{assembly_id}.snippy.fna"

	if skip_if_done(outfile, force, quiet=True):
		return "skipped"

	tmp_ctg = snippy_dir / f"{assembly_id}.tmp.fna"

	# Decompress or copy
	try:
		if str(assembly_path).endswith(".gz"):
			with gzip.open(assembly_path, "rb") as f_in, open(tmp_ctg, "wb") as f_out:
				shutil.copyfileobj(f_in, f_out)
		else:
			shutil.copy(assembly_path, tmp_ctg)

		cmd = (
			f"snippy --outdir {snippy_dir}/{assembly_id} "
			f"--ref {ref_genome} --ctgs {tmp_ctg} --force"
		)
		run_container("snippy", cmd, quiet=True)

		# look for possible outputs
		snps_path = snippy_dir / assembly_id / "snps.consensus.subs.fa"
		if snps_path.exists():
			shutil.move(snps_path, outfile)
			return "ok"
		return "warn"

	except Exception as e:
		return f"error: {e}"

	finally:
		if Path(snippy_dir / assembly_id).exists():
			shutil.rmtree(snippy_dir / assembly_id, ignore_errors=True)
		if tmp_ctg.exists():
			tmp_ctg.unlink(missing_ok=True)


def run_snippy_all(df, snippy_dir: Path, ref_genome: Path, threads: int, force=False):
	ensure_dir(snippy_dir)

	n_total = len(df)
	results = defaultdict(int)  # dynamically count any status

	print(f"[INFO] Running Snippy on {n_total} assemblies (threads={threads})...")

	def process_row(row):
		return row["assembly"], run_snippy_single(row["assembly"], Path(row["assembly_path"]), snippy_dir, ref_genome, force=force)
	
	def print_progress_line():
		# print all status counts on one line
		total_done = sum(results.values())
		status_parts = [f"{k.upper()}: {v}" for k, v in results.items() if k != "error"]
		line = f"[SNIPPY] {' | '.join(status_parts)} (done {total_done}/{n_total})"
		sys.stdout.write(f"\r{line}")
		sys.stdout.flush()

	if threads == 1:
		# sequential
		for _, row in df.iterrows():
			asm, status = process_row(row)
			results[status.split(":")[0]] += 1
			if key == "error":
				results[key] += 1
				print(f"\n[ERROR] {asm} failed")
			else:
				results[key] += 1
				print_progress_line()
	else:
		# multithreaded
		with ThreadPoolExecutor(max_workers=threads) as executor:
			futures = {executor.submit(process_row, row): row["assembly"] for _, row in df.iterrows()}
			for future in as_completed(futures):
				asm = futures[future]
				try:
					_, status = future.result()
					key = status.split(":")[0]
					if key == "error":
						results[key] += 1
						print(f"\n[ERROR] {asm} failed")
					else:
						results[key] += 1
						print_progress_line()
				except Exception as e:
					results["error"] += 1
					print(f"[ERROR] {asm} failed: {e} (done {sum(results.values())}/{n_total})")

	# final summary
	print("\n[SUMMARY] Snippy finished:")
	for status, count in results.items():
		print(f"  {status.upper()}: {count}")
	print(f"  TOTAL: {sum(results.values())}/{n_total}")

# ============================================================
# Step 5: Concatenate Snippy results
# ============================================================

def concat_snippy(df, snippy_dir: Path, snippy_out: str, force=False):
	"""
	Concatenate the longest contig from each per-assembly Snippy file into a single FASTA.
	The output FASTA format is simplified: one line per sequence.
	"""
	if skip_if_done(snippy_out, force):
		return
	print(f"[INFO] concatenating snippy results")
	snippy_files = []
	for _, row in df.iterrows():
		assembly_id = row["assembly"]
		snippy_file = snippy_dir / f"{assembly_id}.snippy.fna"
		if snippy_file.exists():
			snippy_files.append((assembly_id, snippy_file))
		else:
			print(f"[WARN] Snippy output not found for {assembly_id}, skipping")

	if not snippy_files:
		print("[WARN] No Snippy files found to concatenate.")
		return
	
	print(f"[INFO] {len(snippy_files)} are to be concatenated")

	try:
		with open(snippy_out, "w") as out_f:
			for assembly_id, f in snippy_files:
				# Find the longest contig
				contigs = list(SeqIO.parse(f, "fasta"))
				if not contigs:
					print(f"[WARN] No contigs in {f}, skipping")
					continue
				longest_contig = max(contigs, key=lambda x: len(x.seq))
				# Write in single-line format
				seq_line = str(longest_contig.seq)
				out_f.write(f">{assembly_id}\n{seq_line}\n")
	except Exception as e:
		print(f"[ERROR] Failed to write concatenated file: {e}")
		return
	print(f"[OK] Concatenated longest contigs → {snippy_out} ({len(snippy_files)} assemblies)")
	return 
	# Only remove individual files after successful concatenation
	for _, f in snippy_files:
		try:
			f.unlink()
		except Exception as e:
			print(f"[WARN] Could not remove {f}: {e}")

	





# ============================================================
# Step 6: Remove duplicates using seqkit
# ============================================================

def remove_duplicates(snippy_out: Path, projectID: str, threads: int, force=False):
	outfile = snippy_out.with_suffix(".nodup.fna")
	if skip_if_done(outfile, force):
		return outfile
	cmd = f"seqkit rmdup -s {snippy_out} -D {snippy_out.parent}/duplicates.txt -o {outfile} -j {threads}"
	run_container("seqkit", cmd)
	print(f"[OK] Duplicates removed → {outfile}")
	return outfile


# ============================================================
# Step 7: Run Gubbins
# ============================================================

def run_gubbins(snippy_nodup: Path, projectID: str, threads: int, force=False):
	outdir = snippy_nodup.parent
	snps = outdir / f"{projectID}.snps.fasta"
	tree = outdir / f"{projectID}.gubbins.nwk"
	
	if skip_if_done(snps, force) and skip_if_done(tree, force):
		return snps, tree

	cmd = f"run_gubbins.py --model-fitter raxmlng --tree-builder fasttree --threads {threads} --iterations 10 {snippy_nodup}"
	try:
		result = run_container("gubbins", cmd)
	except Exception as e:
		print(f"[ERROR] Gubbins execution failed: {e}")
		raise RuntimeError(f"Gubbins failed for {snippy_nodup}")

	# Check if expected output files exist
	missing = []
	for f in [snps, tree]:
		if not f.exists():
			missing.append(f.name)
	if missing:
		raise RuntimeError(f"[ERROR] Gubbins did not produce expected output files: {', '.join(missing)}")

	print(f"[OK] Gubbins finished → {snps}, {tree}")
	return snps, tree



# ============================================================
# Step 8–11: Tree operations (R-based)
# ============================================================

def run_R(script, args_line, label, container="r_container"):
	print(f"[INFO] Running R step: {label}")
	run_container(container, f"Rscript {script} {args_line}")


def root_tree(project_dir: Path, tree: Path, assemblies: Path, threads: int, params, label: str, force=False):
	outfile = tree.with_suffix(f".rooted_{label}.nwk")
	if skip_if_done(outfile, force):
		return outfile
	cmd = (
		f"--project_dir {project_dir} --tree {tree} --assemblies {assemblies} "
		f"--root_method {params['root_method']} --root_topn {params['root_topn']} "
		f"--threads {threads} --output {outfile}"
	)
	run_R(f"{project_dir}/bin/root_tree.R", cmd, f"root_{label}", project_dir)
	return outfile


def date_tree(project_dir: Path, rooted_tree: Path, snps: Path, assemblies: Path, threads: int, params, label: str, force=False):
	outfile = rooted_tree.with_suffix(f".dated_{label}.rds")
	if skip_if_done(outfile, force):
		return outfile
	cmd = (
		f"--project_dir {project_dir} --trees {rooted_tree} --snps {snps} "
		f"--assemblies {assemblies} --threads {threads} "
		f"--branch_dimension snp_per_genome --clock {params['clock']} --reroot false"
	)
	run_R(f"{project_dir}/bin/date_tree.R", cmd, f"date_{label}", project_dir)
	return outfile


def choose_tree(project_dir: Path, dated_shrink: Path, dated_prune: Path, force=False):
	outfile = project_dir / "best_tree_selected.rds"
	if skip_if_done(outfile, force):
		return outfile
	cmd = f"--trees_shrink {dated_shrink} --trees_prune {dated_prune}"
	run_R(f"{project_dir}/bin/choose_dated_tree.R", cmd, "choose_best", project_dir)
	return outfile


def add_duplicates(project_dir: Path, launch_dir: Path, tree: Path, duplicates: Path, force=False):
	outfile = project_dir / "final_tree_with_duplicates.nwk"
	if skip_if_done(outfile, force):
		return outfile
	cmd = (f"--project_dir {project_dir} --launch_dir {launch_dir} --tree {tree} --duplicates {duplicates}")
	run_R(f"{project_dir}/bin/add_duplicates.R", cmd, "add_dups", project_dir)
	return outfile


# ============================================================
# Main
# ============================================================


def setup_workdir(projectID: str) -> Path:
	workdir = Path(f"treebuilder/{projectID}_work")
	ensure_dir(workdir)
	return workdir

def load_and_validate_assemblies(assemblies_file: Path, workdir: Path, force=False) -> pd.DataFrame:
	df = read_assemblies(assemblies_file)
	validate_assemblies(df, workdir, force)
	return df

def prepare_reference(df: pd.DataFrame, workdir: Path, projectID: str, force=False) -> Path:
	ref_genome = choose_reference_genome(df, workdir, projectID, force)
	return ref_genome

def run_snippy_pipeline(df: pd.DataFrame, snippy_dir: Path, ref_genome: Path, projectID: str, threads: int, force=False) -> Path:
	snippy_out = snippy_dir.parent / f"{projectID}.snippy_out.fna"
	ensure_dir(snippy_dir)

	if skip_if_done(snippy_out, force):
		print(f"[INFO] Snippy output already exists → {snippy_out}, skipping Snippy step.")
		return snippy_out

	run_snippy_all(df, snippy_dir, ref_genome, threads, force)
	concat_snippy(df, snippy_dir, snippy_out, force)
	return snippy_out

def run_downstream_pipeline(snippy_out: Path, projectID: str, threads: int, workdir: Path, assemblies_file: Path, args):
	params = {
		"root_method": args.root_method,
		"root_topn": args.root_topn,
		"clock": args.clock
	}
	
	snippy_nodup = remove_duplicates(snippy_out, projectID, threads)
	snps, gubbins_tree = run_gubbins(snippy_nodup, projectID, threads)

	tree_shrink = workdir / f"{projectID}.treeshrink.nwk"
	tree_prune = workdir / f"{projectID}.treepruner.nwk"

	rooted_shrink = root_tree(workdir, tree_shrink, assemblies_file, threads, params, "shrink")
	rooted_prune = root_tree(workdir, tree_prune, assemblies_file, threads, params, "prune")

	dated_shrink = date_tree(workdir, rooted_shrink, snps, assemblies_file, threads, params, "shrink")
	dated_prune = date_tree(workdir, rooted_prune, snps, assemblies_file, threads, params, "prune")

	best_tree = choose_tree(workdir, dated_shrink, dated_prune)
	final_tree = add_duplicates(workdir, Path.cwd(), best_tree, workdir / "duplicates.txt")

	return final_tree


def main():
	args = get_args()

	# Step 1-3: Setup & prepare assemblies
	workdir = setup_workdir(args.projectID)
	df = load_and_validate_assemblies(Path(args.assemblies), workdir, args.force)
	ref_genome = prepare_reference(df, workdir, args.projectID, args.force)

	# Step 4-5: Snippy
	snippy_dir = workdir / "snippy"
	snippy_out = run_snippy_pipeline(df, snippy_dir, ref_genome, args.projectID, args.threads, args.force)

	# Step 6-11: Downstream analysis
	final_tree = run_downstream_pipeline(snippy_out, args.projectID, args.threads, workdir, Path(args.assemblies), args)

	print(f"[DONE] Pipeline complete → {final_tree}")


if __name__ == "__main__":
	main()
