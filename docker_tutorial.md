# BMEB Bootcamp 2025 Using and creating Containers

Tutorial loosely based on [https://docs.sylabs.io/guides/2.6/user-guide/quick\_start.html](https://docs.sylabs.io/guides/2.6/user-guide/quick_start.html)

Containers are great to familiarize yourself with as a bioinformatician, as they can solve many issues with reproducibility. This tutorial only goes over the absolute basics so you are familiar with using containers if you haven't before.

## Environments - Conda

First let's allocate some resources for our session and load singularity

```
salloc --partition=instruction -N 1 -n 1 -t 03:00:00 
ssh ${SLURM_NODELIST}

module load singularity-ce/singularity-ce.4.1.4
```

## Pulling and building an existing image with Singularity

```
singularity pull shub://vsoch/hello-world 

singularity build hello-world.simg shub://vsoch/hello-world
```

Now we can interact with our container via an interactive session with`singularity shell`

```
singularity shell hello-world.simg

# try running the following commands 
whoami
id

# to leave your shell you can just use 
exit
```

We can also run container without first entering into it with shell

```
singularity run hello-world.simg
```

You can also use run with containers directly from shub and docker

```
singularity run shub://GodloveD/lolcow
```

Here are a few more useful commands you can use to check out your containers

```
singularity inspect hello-world.simg
```

If we want to just run a single command within our container we can use `singularity exec`

```
singularity exec hello-world.simg ls
```

### tricky things with containers

One of the most common issues when working with containers is not mounting your files correctly. This means that you enter into your container and can't find the datafiles you wanted to work with. Singularity automatically will mount your execution directory, home directory, and /tmp, but you may want to mount other locations like this

```
# first let's create or select a file that we want to mount
mkdir ~/testData

nano ~/testData/data.txt
# write some text into the file, it doesn't matter what 

# now we can mount that directory when we execute a command in our container 

singularity exec --bind ~/testData:/mnt hello-world.simg cat /mnt/data.txt
```

That exec command should show you the contents of data.txt

### finding containers

There are a number of container repositories that are compatible across different container platforms. Often when you need a container there is an existing one you can use online. Here are some places you can find containers
[https://biocontainers.pro/registry](https://biocontainers.pro/registry)
[https://quay.io/](https://quay.io/)
[https://hub.docker.com/](https://hub.docker.com/)

Now go find a container that has some software you are familiar with, pull the image, and see if you can view the software documentation with a `singularity exec` command

It should look something like this

```
singularity exec image.simg softwareName --help  
```

# 🧬 Bioinformatics Workflow Example with Singularity

Now that you know how to pull and run containers, let’s run a **mini bioinformatics pipeline** using BioContainers.

We’ll go through three steps:

1. **FastQC** – quality control of reads
2. **BWA** – aligning reads to a reference genome
3. **SAMtools** – converting, sorting, and generating alignment stats

This will give you a taste of how containers are used in real research workflows.

## Step 1: Quality Control with FastQC

Pull the container:

```
singularity pull docker://biocontainers/fastqc
```

Run FastQC on a FASTQ file:

```
singularity exec fastqc_latest.sif fastqc sample.fastq
```

**Expected output:**

* `sample_fastqc.html` (interactive QC report)
* `sample_fastqc.zip` (raw QC metrics)

👉 Open the HTML file to explore base qualities, GC content, and read length distribution.

## Step 2: Alignment with BWA

Pull the container:

```
singularity pull docker://biocontainers/bwa
```

Run alignment:

```
singularity exec bwa_latest.sif bwa mem reference.fasta sample.fastq > aligned.sam
```

**Checkpoint:**

* Inspect first lines:

```
head aligned.sam
```

Notice the header lines (`@SQ`, `@PG`) and read alignments.

## Step 3: Alignment Stats with SAMtools

Pull the container:

```
singularity pull docker://biocontainers/samtools
```

Convert and sort:

```
singularity exec samtools_latest.sif samtools view -bS aligned.sam > aligned.bam
singularity exec samtools_latest.sif samtools sort aligned.bam -o aligned.sorted.bam
singularity exec samtools_latest.sif samtools index aligned.sorted.bam
```

Get alignment stats:

```
singularity exec samtools_latest.sif samtools flagstat aligned.sorted.bam
```

**Expected output:**

* Summary of how many reads mapped/unmapped.

## Wrap-Up

In this tutorial you learned how to:

* Pull existing containers from online registries
* Run bioinformatics tools inside Singularity on HPC
* Build a simple reproducible workflow:
  **FastQC → BWA → SAMtools**

Next steps:

* Try swapping in your own FASTQ + reference genome
* Explore other tools in BioContainers (STAR, bedtools, etc.)
* Think about how to link multiple containerized tools into full pipelines (Snakemake, Nextflow)

[Docker](https://docker-curriculum.com/)
[Create your own container image with docker](https://chtc.cs.wisc.edu/uw-research-computing/docker-build)
