SQANTI3 can be run in any machine that can install python, conda and the necessary dependencies. However, the time complexity and memory requirements can be very different depending on the dataset analyzed.

To test the evolution of time and memory complexity as the dataset grows, we used one medium size transcriptome with 156k transcripts and a big transciptome with 3.56 million transcripts.

# Time complexity (no parallelization)

The necessary time to run SQANTI3 (QC, filter or rescue) increases linearly with the number of transcripts given in the input file (either .gtf or fastq). Therefore, SQANTI3 QC has a time complexity that differs in different stages:

* As the number of transcripts in the gtf or fastq input file increases, the time needed increases linearly (time complexity is **O(n)**), 
* The longest step is usually the mapping of the short reads, which can take most of computation time if the transcriptome is small in comparison to the short reads files. To speed up the process, we recommend mapping the reads beforehand, and giving the mapped bam as input.

# Memory complexity (no parallelization)

SQANTI3 can be used in most desktops and laptops for small and medium transcriptomes: For example, a 156k transcripts require about 5.8 GB of RAM, therefore it can be run on a desktop or laptop with 8 GB of RAM. Regardless, it may take a considerable amount of time, as this example transcriptome took 28 minutes to run in an HPC environment.

In any case, some of the parameters used to run SQANTI3 require and increase in the amount of RAM, such as providing a short reads to map.


# Paralellization requirements

Using the -n option in chunks allows SQANTI3 QC to work in parallel and accelerate the analysis. However, by doing so, memory requirements increase in a considerable manner.

For example, using the 156k transcripts dataset for testing parallelization process in a HPC, the memory cost is linear (memory complexity is **O(n)**, n being how many chunks were used), while real time needed to process the dataset was reduced up to a 25% of the time needed to process without parallelization.

* This dataset took 1678 seconds (about 28 minutes) to process without parallelization, and required 5.87 GB of RAM.
* Using -n 2 reduced the time to 624 seconds (about 10 minutes), and increased memory to 8.12 GB. 
* Using -n 7 further reduced the time to 525s (about 8 and half minutes).
* Using -n 14 (the optimal number for this dataset), took 322 seconds and required 48 GB of RAM.
* From this moment onwards, adding more chunks improves the run time marginally, while RAM usage still grows linearly. In cases where memory is not a constraint and time is limited, it may be useful to use as many chunks as your memory can hold.

Time is greatly reduced just by using parallelization.

![image](https://github.com/user-attachments/assets/91b4a1e2-84bb-4db9-a885-f5b7360f6a8d)

However, memory increase must be taken into account if memory might be a limitation in analysis. The linear increase may be a problem in some use cases.

![image](https://github.com/user-attachments/assets/f2cc8b02-b8e3-4370-825e-d4c81434c227)


For the big (3.56 million transcripts) dataset, time and memory greatly increased:

* This dataset took ~22 hours to be processed without parallelization, requiring 23.89 GB of RAM. With this RAM usage, transcriptomes as big as this one might not be able to be processed in desktop computers or laptops.
* Using -n 4 or -n 5, reduced the run time to approximately 10 hours, but the memory increased up to 100 GB of RAM.
* Using -n 10 required 200 GB of RAM.

# Outliers

When testing the use of different number of chunks, some outliers were found in the processing time of the dataset. In the medium dataset, running SQANTI3 QC using 5, 6 or 15 chunks took longer than using 4 or 14 chunks, while memory usage increased. These results were consistent across several repetitions. While the general trend does not change (a higher number of chunks reduces the overall time required, but increases memory usage in a linear manner), this suggests that some trial-and-error might be needed to find the optimal number of chunks for each dataset.

# Conclusions and recommendations

Increasing the number of parallelization chunks can speed SQANTI3 QC to take as low as 25% of the non-parallelized time. However, the linear increase in memory requirements implies that this option should be used with care in situations in which memory might be a constraint. If SQANTI3 must be used in a desktop, laptop or server with limited specifications, parallelization should be planned in small datasets, which may be fast enough even without parallelization, or medium-size, with a reduced number of cores. Medium or big datasets should be analyzed with 32 GB of RAM or more, but parallelization can be risky, and out of memory errors may be common.

Big datasets, such as the one mentioned above, are recommended to be executed in HPC systems if available, as RAM requirements and time needed to analyze them increase considerably. Using 5 or more cores would require at least 100 GB of RAM.

Regardless of the situation, similar transcriptomes (in number of transcripts) should have a similar performance in time and memory when analyzed with the same parameters. Therefore, once you are confident the transcriptome is analyzed with a good, stable number of cores, you can consider that other similar datasets should have similar requirements. This can be used to estimate how much RAM is needed and how many chunks should be used.