Using Snakemake
===

This is adapted from the Bioinformatics community of practice [HTC and Pipelines](https://github.com/kipkurui/hpc_and_cloud_computing_pipelines/blob/master/snakemake_tutorial.md) course. 

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a command-line oriented workflow management system. In some ways it is conceptually similar to the Linux [make](https://www.tutorialspoint.com/unix_commands/make.htm) command in that it works using filename patterns and reasons about how to produce outputs given a collection of inputs and rules.

## Installing Snakemake

With the Bioconda channel installed, Snakemake can be installed using Bioconda with the `conda install snakemake` command.

## A first Snakemake example

Snakemake expects to be run from a working directory containing a file named `Snakefile` which contains [rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html) describing how to create an output from an input. For this example imagine an input file in the path `data/text.txt` containing the text `gattaca`. Then imagine a `Snakefile` with the content:

```
rule complement:
  input:
    "data/text.txt"
  output:
    "complement_data/text.txt"
  shell:
    "cat {input}|tr atcg tagc > {output}"
```

The directory layout should now look like:

```
data/text.txt
Snakefile
```

This `Snakefile` contains a single rule that explains how to create a file called `complement_data/text.txt` from the input `data/text.txt`. It is very simple and very static, but allows for the action of the rule (the *shell*) bit to be tested. To visualise what the rule does, run `snakemake` with the `-np`:

```
$ snakemake -np 
   input: data/text.txt
    output: complement_data/text.txt
    jobid: 0

cat data/text.txt|tr atcg tagc > complement_data/text.txt
Job counts:
        count   jobs
        1       complement
        1
```

In the output you can see that `snakemake` plans to run the *complement* rule and that it will run the shell command:

```
cat data/text.txt|tr atcg tagc > complement_data/text.txt
```

which runs **tr** which in turn takes the output from **cat** and substitutes `t` for`a`, `a` for `t` and so forth, within input coming from `data/text.txt` and output going to `complement_data/text.txt`. Notice how these file names (from the `input:` and `output:` clauses in the rule) have taken the places of the `{input}` and `{output}` places in the `shell:` command. After examining the working of the rule, we can run it:

```
$snakemake

Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       complement
        1

rule complement:
    input: data/text.txt
    output: complement_data/text.txt
    jobid: 0

Finished job 0.
1 of 1 steps (100%) done
```

and now the directory structure looks like:

```
data/text.txt
complement_data/text.txt
Snakefile
```
with `complement_data/text.txt` containing `ctatggt`, i.e. the complementary DNA bases to those contained in `data/text.txt`. Notice that `snakemake` created the output directory so that the shell command would have a location to write data to.

## Parameterising rules with wildcards

The rule explained above is a static rule. It always does the same action. That is not so useful for real world applications as it means you need to change the `Snakefile` to change what happens. The static rule was not useless though, it illustrated that the shell command works and was thus a useful incremental step in development.

To make the rule re-useable we need to give it parameters. Here is the updated `Snakefile`:

```
rule complement:
  input:
    "data/{file}"
  output:
    "complement_data/{file}"
  shell:
    "cat {input}|tr atcg tagc > {output}"
```

The filename parts of the `input:` and `output:` clause are now replaced with the wildcard (or placeholder) `{file}`. This also means that the rule cannot be run without parameters, as `snakefile -np` illustrates:

```
$ snakemake -np
WorkflowError:
Target rules may not contain wildcards. Please specify concrete files or a rule without wildcards.
```

To use this `Snakefile` we need to give a concrete form to the output. Before proceeding, please *delete* the old output (`rm -r complement_data`) to ensure that `snakemake` has some work to do. Now that the working directory is clean, run `snakemake complement_data/test.txt`:

```
$ snakemake complement_data/test.txt
iMssingInputException in line 1 of /home/user/Workflows_and_package_management/Snakefile:
Missing input files for rule complement:
data/test.txt
```

Oops! Another error. Here the problem is that we asked for a file named `complement_data/test.txt` so the `{file}` wildcard became `test.txt` and there is no file in `data/` named `test.txt` to use as input. The correct command is `snakemake complement_data/text.txt`:

```
$ snakemake complement_data/text.txt
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       complement
        1

rule complement:
    input: data/text.txt
    output: complement_data/text.txt
    jobid: 0
    wildcards: file=text.txt

Finished job 0.
1 of 1 steps (100%) done
```
With the correct wildcard `text.txt` the rule can now run.

### Declarative maps and the DAG

Why specify the *output* file rather than the input file though? Like its predecessor `make`, `snakemake` *reasons* backwards from the required output. Many programming languages are *imperative*, that it programs in languages like R and Python involve specifying what actions the computer needs to take. A `Snakefile` is *declarative*. It is like a map showing which roads head in which direction. The inputs to `snakemake` are a destination (or result), the rules of what moves are permissible, and of course the current state of the work directory. Given this information `snakemake` tries to plot a route to the desired output.

This is a bit like being at being at [KEMRI-WT](https://www.kemri-wellcome.org/) and wantin to travel. It is not enough to state that you are at ILRI and want to move. You need to specify a destination (for example Malindi) and then a good enough map will give you directions to get there.


## Expanding beyond a single rule

Workflows naturally are composed of more than single step or rule. In the output from `snakemake` there is a reference to a DAG. That is a [Directed Acyclic Graph](https://en.wikipedia.org/wiki/Directed_acyclic_graph). Leaving aside the computer science of DAGs, they are the form that the *rule map* takes. Each rule is a node in the graph and connections between rules (which will be shown soon) are the edges.

This model corresponds to how one might sketch a workflow on a whiteboard. Each action (*rule*) transforms data from one form to the next, and the flow of data between rules would be the arrows that you draw between actions. In fact a classic model of bioinformatics suggests that it is all [transforms and filters](https://academic.oup.com/bioinformatics/article/30/17/i601/201260). In `snakemake` what binds rules together is common filenames, as can be seen in the following `Snakefile`:

```
rule complement:
  input:
    "data/{file}"
  output:
    "complement_data/{file}"
  shell:
    "cat {input}|tr atcg tagc > {output}"

rule reverse:
  input:
    "complement_data/{file}"
  output:
    "reversed_data/{file}"
  shell:
    "cat {input} |rev > {output}"
```

The previously described `complement` rule appears again, and now there is a second rule, `reverse` that uses the *rev* command to reverse the a line in a file. Together the two rules implement the [reverse complement](http://genewarrior.com/docs/exp_revcomp.jsp) operation so common in bioinformatics. As mentioned previously the filenames bind the rules together, in this case the `complement_data/{file}` pattern common to the `output:` of `complement` and the `input:` of `reverse`. As before, remove the previous output (`rm -r complement_data`). The output we now desire is `reversed_data/text.txt` and to get the work place run `snakemake -np reversed_data/text.txt`:

```
$ snakemake -np reversed_data/text.txt
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       complement
        1       reverse
        2

rule complement:
    input: data/text.txt
    output: complement_data/text.txt
    jobid: 1
    wildcards: file=text.txt

Finished job 1.
1 of 2 steps (50%) done

rule reverse:
    input: complement_data/text.txt
    output: reversed_data/text.txt
    jobid: 0
    wildcards: file=text.txt

Finished job 0.
2 of 2 steps (100%) done
```

Which shows us that 1 `complement` job is needed and 1 `reverse` job, for a total of 2 jobs. If you have the `gnuplot` and `ImageMagick` software installed you can also get a graphical overview of the work plan with the command `snakemake --dag reversed_data/text.txt | dot | display`.

After running the command `snakemake reversed_data/text.txt` we the working directory now contains:

```
complement_data/text.txt
data/text.txt
reversed_data/text.txt
Snakefile
```
with `reversed_data/text.txt` containing `tggtatc`, i.e. the reverse complement of the starting sequence `gatacca`.

This concludes the reverse complement example. You are urged to have a look at the official [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) and to continue to the next example which illustrates a bit of *software dependency management*, a crucial aspect of reproducible and reliable scientific workflows.

## Working with dependencies

The reverse complement example used tools that are installed with (almost) every Linux installation. This is not the case for most of our research workflows. For the next example we will need two data files in a directory called 'data/'. We can download them from Zenodo as follows:

```
wget -O data/A.fastq.gz http://bit.ly/sm_r1_fq
wget -O data/B.fastq.gz http://bit.ly/sm_r2_fq
```

and then the `Snakefile` is as follows:

```
rule fastqc:
  input:
    "data/{sample}.fastq.gz"
  output:
    "{out_dir}/{sample}_fastqc.zip",
    "{out_dir}/{sample}_fastqc.html"
  conda:
    "fastqc.yml"
  shell:
    "fastqc --outdir {wildcards.out_dir} {input}"
```

This `Snakefile`, designed to run the [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) quality reporting program, has two features not previously shown. Firstly, it uses a wildcard for the output directory, '{out_dir}'. This is used in the `shell:` clause where it is referred to as `{wildcards.out_dir}`. The `wildcards` prefix allows wildcard patterns to be specified and used as variables in `shell:` expression.

Secondly it has a `conda:` clause. This references to a file that describes the [conda](https://bioconda.github.io/) environment that needs to be created to provide the software dependencies this rule needs. That `fastqc.yml` file was created using the `conda env export` command and contains:

```
name: fastqc
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - fastqc=0.11.7=5
  - libgcc-ng=7.2.0=hdf63c60_3
  - perl=5.26.2=h470a237_0
  - openjdk=8.0.152=h46b5887_1
```

If not activated the `conda:` clause is ignored. It requires the `--use-conda` to become active. The `rule` can be tested by requesting some output:

```
$ snakemake -np --use-conda qc_reports/A_fastqc.html
Building DAG of jobs...
Conda environment fastqc.yml will be created.
Job counts:
	count	jobs
	1	fastqc
	1

rule fastqc:
    input: data/A.fastq.gz
    output: qc_reports/A_fastqc.zip, qc_reports/A_fastqc.html
    jobid: 0
    wildcards: sample=A, out_dir=qc_reports

fastqc --outdir qc_reports data/A.fastq.gz
Job counts:
	count	jobs
	1	fastqc
	1
```

As before we are presented a goal to `snakemake` and it responded with a description of the command line that would be run. In addition, it stated that a conda environment would be created. In practice conda environments are created by `snakemake` in subdirectories of the `.snakemake/` directory in the current directory. Running `snakemake --use-conda qc_reports/A_fastqc.html` results in:

```
$ snakemake --use-conda qc_reports/A_fastqc.html
Creating conda environment fastqc.yml...
Environment for fastqc.yml created (location: .snakemake/conda/740e2bc8)
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       fastqc
        1

rule fastqc:
    input: data/A.fastq.gz
    output: qc_reports/A_fastqc.zip, qc_reports/A_fastqc.html
    jobid: 0
    wildcards: out_dir=qc_reports, sample=A

Started analysis of A.fastq.gz
Approx 5% complete for A.fastq.gz
Approx 10% complete for A.fastq.gz
Approx 15% complete for A.fastq.gz
Approx 20% complete for A.fastq.gz
Approx 25% complete for A.fastq.gz
Approx 30% complete for A.fastq.gz
Approx 35% complete for A.fastq.gz
Approx 40% complete for A.fastq.gz
Approx 45% complete for A.fastq.gz
Approx 50% complete for A.fastq.gz
Approx 55% complete for A.fastq.gz
Approx 60% complete for A.fastq.gz
Approx 65% complete for A.fastq.gz
Approx 70% complete for A.fastq.gz
Approx 75% complete for A.fastq.gz
Approx 80% complete for A.fastq.gz
Approx 85% complete for A.fastq.gz
Approx 90% complete for A.fastq.gz
Approx 95% complete for A.fastq.gz
Approx 100% complete for A.fastq.gz
Analysis complete for A.fastq.gz
Finished job 0.
1 of 1 steps (100%) done
```

As the log file shows the dependencies specified in `fastqc.yml` are downloaded and installed into a conda environment and this conda environment is activated to provide `fastqc` and allow the analysis to run. After the run the working directory contains:

```
data/A.fastq.gz
data/B.fastq.gz
qc_reports/A_fastqc.html
qc_reports/A_fastqc.zip
fastqc.yml
Snakefile
```
Note that `fastqc` generated two outputs, which is standard for this tool and is described in the rule.

In this example there are two input files, `A.fastq.gz` and `B.fastq.gz`. These can both be processed by the `fastqc` rule, but for that to happen a report needs to be specified as an output for each tool. To prepare for this, delete the `qc_reports` directory (`rm -r qc_reports/`) and then run `snakemake --use-conda qc_reports/A_fastqc.html qc_reports/B_fastqc.html`. This time the conda environment will not have to be created from scratch as it already exists.

```
$ snakemake --use-conda qc/A_fastqc.html qc/B_fastqc.html
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        2       fastqc
        2

rule fastqc:
    input: data/A.fastq.gz
    output: qc_reports/A_fastqc.zip, qc_reports/A_fastqc.html
    jobid: 0
    wildcards: out_dir=qc_reports, sample=A

Started analysis of A.fastq.gz
Approx 5% complete for A.fastq.gz
Approx 10% complete for A.fastq.gz
Approx 15% complete for A.fastq.gz
...
...
...
Approx 85% complete for A.fastq.gz
Approx 90% complete for A.fastq.gz
Approx 95% complete for A.fastq.gz
Approx 100% complete for A.fastq.gz
Analysis complete for A.fastq.gz
Finished job 0.
1 of 2 steps (50%) done

rule fastqc:
    input: data/B.fastq.gz
    output: qc_reports/B_fastqc.zip, qc_reports/B_fastqc.html
    jobid: 1
    wildcards: out_dir=qc_reports, sample=B

Started analysis of B.fastq.gz
Approx 5% complete for B.fastq.gz
Approx 10% complete for B.fastq.gz
Approx 15% complete for B.fastq.gz
...
...
Approx 85% complete for B.fastq.gz
Approx 90% complete for B.fastq.gz
Approx 95% complete for B.fastq.gz
Approx 100% complete for B.fastq.gz
Analysis complete for B.fastq.gz
Finished job 1.
2 of 2 steps (100%) done
```

The output is trimmed to focus on the essentials. Given two outputs, the rule caused two jobs to be created. Most computers have more than a single CPU core, so `snakemake` has support for running jobs in parallel, using more than one processing core at a time. Once again delete the previous output (`rm -r qc`) and enable parallel processing with the `-j` option, so the command becomes `snakemake -j 2 --use-conda qc/A_fastqc.html qc/B_fastqc.html`:

```
$ snakemake -j 2 --use-conda qc/A_fastqc.html qc/B_fastqc.html
Provided cores: 2
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	2	fastqc
	2

rule fastqc:
    input: data/B.fastq.gz
    output: qc/B_fastqc.zip, qc/B_fastqc.html
    jobid: 0
    wildcards: out_dir=qc, sample=B
    
rule fastqc:
    input: data/A.fastq.gz
    output: qc/A_fastqc.zip, qc/A_fastqc.html
    jobid: 1
    wildcards: out_dir=qc, sample=A

Started analysis of A.fastq.gz
Started analysis of B.fastq.gz
Approx 5% complete for B.fastq.gz
Approx 5% complete for A.fastq.gz
Approx 10% complete for B.fastq.gz
Approx 10% complete for A.fastq.gz
Approx 15% complete for B.fastq.gz
Approx 15% complete for A.fastq.gz
...
Approx 90% complete for A.fastq.gz
Approx 90% complete for B.fastq.gz
Approx 95% complete for A.fastq.gz
Approx 95% complete for B.fastq.gz
Approx 100% complete for A.fastq.gz
Approx 100% complete for B.fastq.gz
Analysis complete for A.fastq.gz
Analysis complete for B.fastq.gz

Finished job 0.
1 of 2 steps (50%) done
Finished job 1.
2 of 2 steps (100%) done
```

Once again the output is trimmed for clarity. It shows that the processing of `A.fastq.gz` and `B.fastq.gz` proceeded in parallel, that is simultaneously.

All of these examples have shown `snakemake` being run directly on the command line. When using a cluster, however, the `snakemake` command can be included in a shell script submitted to the cluster and the `-j` option tuned according to the CPUs requested from the cluster.

In addition to running on a worker node, `snakemake` can also orchestrate tasks running on a cluster. This is part of a much more sophisticated work mode of `snakemake` where resource requirements for individual stages can be specified along with other cluster-specific tuning. In this case a machine needs to be available to run `snakemake` as a workflow orchestrator. Since `snakemake` will not use a lot of resources, perhaps a small virtual machine can be set aside to host running `snakemake` jobs. The snakemake manual describes useful options for [running on a cluster](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) and [cluster configuration](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration).
