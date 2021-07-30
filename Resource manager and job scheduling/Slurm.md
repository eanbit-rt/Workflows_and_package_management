
This Tutorial is based on Icipe HPC cluster [For more information about icipe cluster](https://hpc01.icipe.org/)

# Slurm user guide

SLURM is the queue manager used on the ICIPE HPC cluster. You are required to use SLURM to submit jobs to the cluster.

All the commands presented in this guide are to be used from the host core.cluster.france-bioinformatique.fr

## SLURM partitions and nodes

The ICIPE HPC cluster is organized into TWO SLURM partitions. Each partition gathers a set of compute nodes that have similar usage.

The default partition used by SLURM (fast) contains a huge number of nodes and is suitable for most jobs.

To view all partitions available on the cluster run :

```bash
sinfo
```

*Please note that you may not have the rights to use all partitions*

To view all available nodes run :

```bash
sinfo -Nl
```


## Submitting a job to the cluster

They are two commands to submit a job to the cluster:

`srun` to run jobs interactively
`sbatch` to submit a batch job

### Submit a job using `srun`

To learn more about the `srun` command, see [the official documentation](https://slurm.schedmd.com/srun.html). You can also use `srun -h` to get the options.

#### Usage

The job will start immediately after you execute the `srun` command. The outputs are returned to the terminal. You have to wait until the job has terminated before starting a new job. This works with ANY command.

Example:

```bash
srun hostname
```

Example if an interaction is needed:

```bash
module load r
srun --mem 20GB --pty R
```
- `--pty`: will keep the interaction possible
- `--mem 20GB`: will allow 20GB of memory to your job instead of the 2GB by default

#### N.B Submit interactive jobs as scripts or specify node to submit interactive job

```
module load python
srun -p partition --mem 20GB --pty python
```

We have allocated 20Gb of memory for running python interactively on the terminal
### Submit a job using `sbatch`

To learn more about the `sbatch` command, see [the official documentation](https://slurm.schedmd.com/sbatch.html), or use `sbatch -h`. 

#### Usage

The job starts when resources are available, and only returns the job id.
The outputs are sent to file(s).
`sbatch` ONLY works with shell scripts. The batch script may be given to `sbatch` through a filename on the command line, or if no filename is specified, `sbatch` will read in a script from standard input.

#### Batch scripts rules

The script can contain `srun` commands. Each `srun` is a job step. The script must start with shebang (#!) followed by the path to the interpreter
```bash
#!/bin/bash
#!/usr/bin/env python
```

The execution parameters can be set:

At runtime in the command `sbatch`

```bash
sbatch --mem=40GB bowtie2.sbatch
```

Or within the shell `bowtie2.sbatch` itself

```bash
#!/bin/bash
#
#SBATCH --mem 40GB
srun bowtie2 -x hg19 -1 sample_R1.fq.gz -2 sample_R2.fq.gz -S sample_hg19.sam
```

```bash
sbatch bowtie2.sbatch
```


The scripts can contain slurm options just after the shebang but before the script commands → #SBATCH

**Note** that the syntax `#SBATCH ` is important and doesn't contain any `!` (as in the Shebang)

**Advice:** We recommend setting as many parameters as you can in the script to keep track of your execution parameters for future submission.

### Execution parameters

These parameters are common to the commands `srun` and `sbatch`.

#### Parameters for log

```bash
#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out  # STDOUT file with the Node name and the Job ID
#SBATCH -e slurm.%N.%j.err  # STDERR file with the Node name and the Job ID
```

#### Parameters to control the job

**`--partition=<partition_names>`, `-p`**

Request a specific partition for the resource allocation. Each partition (queue in SGE) have their own limits: time, memory, nodes ...

Cf: `sinfo` to know which partitions are available.

**`--mem=<size[units]>`**

Specify the real memory required per node. The default unit is `MB`
(Default: 2GB)

The job is killed if it exceeds the limit.

Note that you can use the variable `$SLURM_MEM_PER_NODE` in the command line to synchronize the software settings and the resource allocated.

**`--time=<time>`, `-t`**

Set a limit on the total run time of the job allocation.

Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".

#### Parameters for multithreading

**`--cpus-per-task=<ncpus>`, `--cpus`, `-c`**

Request a number of CPUs (default 1)

Note that you can use the variable `$SLURM_CPUS_PER_TASK` in the command line to avoid mistakes between the resource allocated and the job.

```bash
#!/bin/bash
#
#SBATCH --cpus-per-task=8

srun bowtie2 --threads $SLURM_CPUS_PER_TASK -x hg19 -1 sample_R1.fq.gz -2 sample_R2.fq.gz -S sample_hg19.sam
```


### Full example of the `sbatch` command

#### Random

1. Open a script file with any text editor (but not Word)

For beginners, we suggest to use `nano`, which has restricted functionalities but is quite intuitive.

```bash
nano slurm_random.sh
```

2. Copy/Paste the following script, which is for submitting mafft job to the main hpc01 cluster :


```
#!/usr/bin/env bash
#SBATCH -p qbatch
#SBATCH --nodes=1
#SBATCH -w hpc01
#SBATCH --ntasks=1
​
#loading the module mafft
module load mafft/7.475
​
# Executing the alignment
​
#mafft --auto --reorder --preservecase /mnt/nfs/home/admin/Data/Raw_Data/mtDNA46K.fasta > /mnt/nfs/admin/admin/Data/Derived_Data/mtDNA46Kalign.fasta
#use script for refence only 
​
echo "Successfully aligned admin."
```
Press Ctrl-x to exit nano, then "Y" when nano asks you whether the modified buffer should be saved, then press the "Enter" key to confirm the file name.

3. Check the content of the script

```bash
cat slurm_random.sh
```


4. Submit the job

```bash
sbatch slurm_random.sh
```

5. Check the result

Since this script is running a very basic task, the results should promptly be available.

Check the output files with `ls` and `head`.

Note: these commands cannot be run on the login node since they consume massive compute resources.



### Job information

List a user's current jobs:

```bash
squeue -u <username>
```

List a user's running jobs:

```bash
squeue -u <username> -t RUNNING
```

List a user's pending jobs:

```bash
squeue -u <username> -t PENDING
```

View accounting information for all user's job for the current day :

```bash
sacct --format=JobID,JobName,User,Submit,ReqCPUS,ReqMem,Start,NodeList,State,CPUTime,MaxVMSize%15 -u <username>
```

View accounting information for all user's job for the 2 last days (it worth an alias) :

```bash
sacct -a -S $(date --date='2 days ago' +%Y-%m-%dT%H:%M) --format=JobID,JobName,User%15,Partition,ReqCPUS,ReqMem,State,Start,End,CPUTime,MaxVMSize -u <username>
```

List detailed job information:

```bash
scontrol show -dd jobid=<jobid>
```

### Manage jobs

To cancel/stop a job:

```bash
scancel <jobid>
```

To cancel all jobs for a user:

```bash
scancel -u <username>
```

To cancel all pending jobs for a user:

```bash
scancel -t PENDING -u <username>
```
