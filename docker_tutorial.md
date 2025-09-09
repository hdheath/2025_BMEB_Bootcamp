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

Here’s your tutorial revised so it works correctly with **Quay.io BioContainers** and explicit version tags instead of `latest` (which doesn’t exist). I’ve updated FastQC, BWA, and SAMtools accordingly:

---

Perfect — thanks for clarifying. Since in your workshop/tutorial students will **cd into a working directory after pulling the image**, they’ll run into the same “file not found” issue unless we show them how to either (1) reference the `.sif` with its absolute path, or (2) move it into the working directory. And since their data will be in a separate folder (not in your `/hb/home/hdheath/...`), we should also explicitly teach them how to handle paths and mounts.

Here’s how I’d update your tutorial to make this fool-proof:

---

Got it 👍 — here’s the tutorial rewritten so **students use their own working directories** with **placeholders** (`/path/to/...`) instead of your personal `/hb/home/hdheath/...` paths. That way they can adapt it no matter where they’re running.

---

# 🧬 Bioinformatics Workflow Example with Singularity

Now that you know how to pull and run containers, let’s run a **mini bioinformatics pipeline** using BioContainers.

We’ll go through three steps:

1. **FastQC** – quality control of sequencing reads
2. **BWA** – aligning reads to a reference genome
3. **SAMtools** – converting, sorting, and generating alignment stats

This will give you a taste of how containers are used in real research workflows.

---

## Step 0: Set Up Your Working Directory

First, create a clean working folder and copy the required files into it:

```bash
mkdir -p docker_tutorial
cd docker_tutorial

# Copy reference genome, assembly, and mock reads into your workdir
cp /hb/home/hdheath/bmeb_bootcamp_2025/tardigrade_reference_genome.fa .
cp /hb/home/hdheath/bmeb_bootcamp_2025/tardigrade.fastq .
```

Check that the files are in place:

```bash
ls -lh
# should show:
#   tardigrade_reference_genome.fa
#   tardigrade.fastq
```

---

## Step 1: Quality Control with FastQC

Pull the container (this will create a `.sif` file in your current directory):

```bash
singularity pull docker://quay.io/biocontainers/fastqc:0.11.9--0
```

Run FastQC on the mock reads:

```bash
singularity exec /path/to/fastqc_0.11.9--0.sif fastqc tardigrade.fastq
```

**Expected output:**

* `tardigrade_fastqc.html` (interactive QC report)
* `tardigrade_fastqc.zip` (raw QC metrics)

---

## Step 2: Alignment with BWA

Pull the container:

```bash
singularity pull docker://quay.io/biocontainers/bwa:0.7.17--hed695b0_7
```

Run alignment against the reference genome (`.fna` file):

```bash
singularity exec /path/to/bwa_0.7.17--hed695b0_7.sif \
    bwa mem tardigrade_reference_genome.fa tardigrade.fastq > aligned.sam
```

---

## Step 3: Alignment Stats with SAMtools

Pull the container:

```bash
singularity pull docker://quay.io/biocontainers/samtools:1.9--h10a08f8_12
```

Convert, sort, and index:

```bash
singularity exec samtools_1.9--h10a08f8_12.sif samtools view -bS aligned.sam > aligned.bam
singularity exec samtools_1.9--h10a08f8_12.sif samtools sort aligned.bam -o aligned.sorted.bam
singularity exec /samtools_1.9--h10a08f8_12.sif samtools index aligned.sorted.bam
```

Get alignment stats:

```bash
singularity exec samtools_1.9--h10a08f8_12.sif samtools flagstat aligned.sorted.bam
```

**Expected output:**
A summary of how many reads mapped vs unmapped.

---

## 🔧 Common Pitfalls & Fixes

* **`.sif` file not found** → Always use the absolute path to the `.sif`, e.g. `/path/to/fastqc_0.11.9--0.sif`.
* **Input files missing** → Be sure to `cp` all required files into your working directory first.
* **Files in other storage locations (e.g. `/scratch`)** → Bind them explicitly:

  ```bash
  singularity exec --bind /scratch/data:/mnt /path/to/fastqc_0.11.9--0.sif fastqc /mnt/myreads.fastq
  ```

---

## Wrap-Up

In this tutorial you learned how to:

* Copy files into a **working directory**
* Run **FastQC** on mock reads
* Align reads with **BWA**
* Process alignments with **SAMtools**
* Troubleshoot missing images and data paths

Together these steps form a reproducible workflow:
**FastQC → BWA → SAMtools**

---



