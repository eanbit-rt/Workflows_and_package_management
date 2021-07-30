# Workflows and Package Management

Reproducibility and package management techniques: workflow languages (CWL, Snakemake, and Conda). This course introduces some of the approaches for package management and how to create reproducible workflows or pipelines. 

### Competencies 

This session seeks to impart the following competencies:

1. Knowledge and skills: Bioinformatics tools and their usage.
2. Knowledge and Skills: Command line and scripting based computing skills appropriate to the discipline.

### Learning Outcomes

By the end of this session, and the projects that follow, the learner should be able to:

1. Select the best workflow and package managers based on the task at hand
2. Implement a genomic pipeline in at least one workflow manager
3. Set up a reproducible analysis environment

## Outline

- Introduce the high-level concept of workflows and high throughput data analysis
- Hands-on activities for setting up the packages
- Introduce package management and how we can use conda to increase reproducibility with workflows
- Introduce the theory of workflows: with emphasis on one language (say, snakemake)
- Hands-on activities of developing workflows

## Slides

1. [Using Bioconda to streamline software installation for bioinformatics](https://monashbioinformaticsplatform.github.io/bioconda-tutorial/#/)
2. [Workflows and Pipelines](https://docs.google.com/presentation/d/1AbKftgGsod9dvSTzzgqXc879eQkx-w3oaXXFUtKeNc0/edit?usp=sharing)

## Tutorials

1. [Nextflow and Singularity tutorial](https://galaxyuvri-ea.github.io/nextflow-intro/) by Alfred Ssekagiri.
2. [Docker Tutorial](genomic_workflows_docker_tutorial.md) from Mark Wamalwa
2. [Package Management with conda](miniconda_tutorial.md)
2. [Workflow with Snakemake](snakemake_tutorial.md) will provide a quick introduction then we'll dive deeper using [Reproducible Research](https://reproducibility.sschmeier.com/index.html) tutorial.See [this](https://www.biostars.org/p/335903/) tutorial also
3. [Common Workflow language](nextflow_tutorial.md) tutorial. We will not cover this, but we provide links to useful tutorials for you to explore and learn further. Also see [this](https://www.commonwl.org/user_guide/01-introduction/index.html) and this(https://andrewjesaitis.com/2017/02/common-workflow-language---a-tutorial-on-making-bioinformatics-repeatable/) walkthroughs. 
4. [Resource management on HPC](https://github.com/eanbit-rt/Workflows_and_package_management/blob/59e3233e9cd77bdbbec16c4801e9422217c7cb1e/Resource%20manager%20and%20job%20scheduling/Slurm.md)

## Reading resources

Some resources and articles you can make use in this course:

1. [Awesome pipelines](https://github.com/pditommaso/awesome-pipeline): A curated list of pipelines and workflow languages
2. [Existing Workflow systems](https://github.com/common-workflow-language/common-workflow-language/wiki/Existing-Workflow-systems): Computational Data Analysis Workflow Systems
3. Papers:

    - [A review of bioinformatic pipeline frameworks](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429012/)
    - [Spjuth et al. 2015](https://biologydirect.biomedcentral.com/articles/10.1186/s13062-015-0071-8)
    - [A Review of Scalable Bioinformatics Pipeline](https://link.springer.com/article/10.1007/s41019-017-0047-z)
    - [Bioconda Paper](https://sci-hub.tw/10.1038/s41592-018-0046-7)
    - [July 2019: Scalable Workflows and Reproducible Data Analysisfor Genomics](https://link.springer.com/content/pdf/10.1007%2F978-1-4939-9074-0_24.pdf)
