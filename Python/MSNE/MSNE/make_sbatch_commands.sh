#!/usr/bin/env bash
#
# make_sbatch_commands.sh
#
# This script scans all MSNE_*.sh files in the current directory.
# For each .sh file with NO matching output_<jobname>_clusters.csv,
# it writes an sbatch command to submit_all_msne_jobs.sh.
#
# Usage:
#   1) Make this file executable: chmod +x make_sbatch_commands.sh
#   2) Run ./make_sbatch_commands.sh
#   3) Then submit all jobs at once by: ./submit_all_msne_jobs.sh
#
# Reference for sbatch usage:
#   Slurm Workload Manager Documentation (https://slurm.schedmd.com/sbatch.html)

# Output file containing the sbatch commands
OUTFILE="submit_remaining_msne_jobs.sh"

# Start fresh
echo "#!/usr/bin/env bash" > "$OUTFILE"
echo "# Auto-generated script to submit the remaining MSNE jobs." >> "$OUTFILE"
echo >> "$OUTFILE"

# Collect job scripts to run
jobs_to_run=()
for job in MSNE_*.sh; do
    prefix="${job%.sh}"
    out_clust="output_${prefix}_clusters.csv"
    if [[ ! -f "$out_clust" ]]; then
        jobs_to_run+=("$job")
    fi
done

# Write them all to the output file, one per line
for j in "${jobs_to_run[@]}"; do
    echo "sbatch \"$j\"" >> "$OUTFILE"
done

# Make the output file itself executable
chmod +x "$OUTFILE"

# Optionally notify the user
echo "Created $OUTFILE with ${#jobs_to_run[@]} sbatch commands."
echo "Run './$OUTFILE' to submit all those jobs."
