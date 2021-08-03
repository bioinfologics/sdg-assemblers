# HappySorter pipeline

A pipeline to generate haplotype-specific assemblies.

This pipeline requires datastores for PE and long reads, either Nanopore or PacBio. Generating datastores is covered in the SDG documentation but the following commands will generate a paired-end datastore called 'pe.prseq' and a long read datastore called 'long.loseq';

```
sdg-datastore make -1 pe_reads_R1.fastq -2 pe_reads_R2.fastq -t paired -d 1 -n pe -o pe
sdg-datastore make -L long_reads.fastq -t long -n long -o long
```

In addition, you should know the unique 31-mer coverage peak for the PE reads. For a heterozygous genome this will be the first peak in the 31-mer spectra.

This pipeline has four steps:

* **01-dbg_strider.py** - Create a dbg and extend its haplotype-specific unitigs with short and long reads.

* **02-anchors_rtg.py** - Create anchors and a ReadThreadsGraph.

* **03-make_orders.py** - Use the TotalSorter/HappySorters to create and merge the orders. Includes post-order scaffolding.

* **04-make_outputs.py** - Create output files (assembly consensus, etc).


## 01-dbg_strider.py

```
usage: 01-dbg_strider.py [-h] -o OUTPUT_PREFIX -p PAIRED_DATASTORE -l
                         LONG_DATASTORE [-k K] -u UNIQUE_COVERAGE
                         [-c MIN_COVERAGE] [-b DISK_BATCHES]
                         [--low_node_coverage LOW_NODE_COVERAGE]
                         [--low_bubble_coverage LOW_BUBBLE_COVERAGE]
                         [--lr_min_support LR_MIN_SUPPORT] [--lr_snr LR_SNR]
                         [--lr_max_noise LR_MAX_NOISE]
```

The following arguments are required: -o/--output_prefix, -p/--paired_datastore, -l/--long_datastore, -u/--unique_coverage

## 02-anchors_rtg.py

```
usage: 02-anchors_rtg.py [-h] -o OUTPUT_PREFIX -m MIN_COVERAGE -M MAX_COVERAGE
                         -s MAX_ANCHOR_SIZE
```

The following arguments are required: -o/--output_prefix, -m/--min_coverage, -M/--max_coverage, -s/--max_anchor_size

## 03-make_orders.py

```
usage: 03-make_orders.py [-h] -o OUTPUT_PREFIX
```

The following arguments are required: -o/--output_prefix

## 04-make_outputs.py
