#!/bin/bash

# slurm/disBatch wrapper script written by Nick C. (Simons Foundation)

# Set up a pipe to flow tasks from snakemake to disBatch.
TaskPipe=/tmp/dbTaskPipe_${SLURM_JOBID}
rm -f ${TaskPipe}
mkfifo ${TaskPipe}

# Create a script that uses the pipe to send a cromwell task to
# disBatch.  We will tell cromwell via slurm.config to this script
# (and not squeue) to submit a task.
export SubmitTaskScript=./SubmitViaDB_${SLURM_JOBID}.sh
cat > ${SubmitTaskScript} <<EOF
#!/bin/bash

workdir=\$1
outfile=\$2
errfile=\$3
script=\$4
echo "cd \${workdir} ; bash \${script} > \${outfile} 2> \${errfile}" > ${TaskPipe} 
echo "Submitted batch job \$RANDOM\$RANDOM"
EOF
chmod +x ${SubmitTaskScript}

# cromwell needs a fairly recent java.
module purge
module add modules-new/bash
module add slurm openjdk

# Run disBatch in the back ground, reading tasks from the pipe.
export DISBATCH_ROOT=/cm/shared/sw/pkg/flatiron/disBatch/2.0-beta
${DISBATCH_ROOT}/disBatch -p wdlTest_${SLURM_JOBID} ${TaskPipe} &

# Keep the pipe open.
exec 3>${TaskPipe}

# Run cromwell with config file that routes tasks to disBatch.
java -Dconfig.file=mk-slurm.conf -jar \
/mnt/home/mlek/mlek/bin/cromwell-56.jar run \
MultiSampleReadViz_consolidated.wdl \
-o workflow.options \
-i readviz_inputs_consolidated.json \
&> wdl_${SLURM_JOBID}.log

# Close the pipe
exec 3>&-

# Wait for disBatch to finish up all the tasks.
wait

# Clean up. Probably should put this in an exit function.
rm -f ${TaskPipe} ${SubmitTaskScript}




