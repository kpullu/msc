# A MapReduce approach to Genome Alignment
This MSc dissertation presents two approaches to perform distributed sequence alignment of genomic data based on the MapReduce programming paradigm. 
- **MR-BWA** presents an approach in distributing BWA using MapReduce. [BWA](http://bio-bwa.sourceforge.net/) is an industry standard software used for genomic reads alignment.
- **MR-BWT-FM** presents low level optimizations on suffix array and BWT creation which are used to create a custom FM-Index which in turn is used for distributed genome sequence alignment.
 
Both approaches have been tested on an Apache Hadoop cluster on Amazon EC2 cloud consisting of five nodes of the *r3.4xlarge* instance type. [Apache Hadoop v2.7.1](https://hadoop.apache.org/docs/r2.7.1/), [Python v2.7.6](https://www.python.org/download/releases/2.7.6/) , [BWA v0.7.15](https://sourceforge.net/projects/bio-bwa/files/) and [SAM tools v1.4.1](https://sourceforge.net/projects/samtools/files/samtools/1.4.1/) have been installed on each AWS node. MapReduce jobs were executed using the [Hadoop streaming API v2.6.5](https://hadoop.apache.org/docs/r1.2.1/streaming.html).

## Setting up MR-BWA on Hadoop
It is assumed that an Apache Hadoop Cluster has already been setup
1. Apache Hadoop Cluster configuration is updated to have specific memory assigned:
-- yarn-site.xml
--- ```yarn.log-aggregation-enable=true```
--- ```yarn.nodemanager.pmem-check-enabled=false```
--- ```yarn.nodemanager.vmem-check-enabled=false```
--- ```yarn.nodemanager.resource.memory-mb=12288```
--- ```yarn.scheduler.maximum-allocation-mb=12288```
-- mapred-site.xml
--- ```mapreduce.map.memory.mb=12288```
--- ```mapreduce.reduce.memory.mb=12288```
--- ```mapreduce.map.java.opts=-Xmx12288M -Djava.net.preferIPv4Stack=true -XX:NewRatio=8 -XX:+UseNUMA -XX:+UseParallelGC```
--- ```mapreduce.reduce.java.opts=-Xmx12288M -Djava.net.preferIPv4Stack=true -XX:NewRatio=8 -XX:+UseNUMA -XX:+UseParallelGC```
2. Upload and run the provided ```setup_biotools.sh, hd biotools.sh``` on each node. These scripts will install and setup Python, BWA and SAM tools. Any extra Python dependencies are also installed.
3. Upload the provided ```hadoop-streaming-2.6.5.jar``` on the master node and copy it to the *$HADOOP_HOME* folder.
4. Upload all the provided Python scripts in *mr-bwa* folder to the master node
5. Copy the uploaded Python scripts to *app_mrbwa* folder on the master node. This is needed because Python import modules work with symlink, hence files need to reside in a common folder
6. Create */data/index* folder on master node and download human reference genome in this folder by running the command
```wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz```
7. Uncompress the downloaded human reference genome and index it by running the commands
```gunzip hg38.fa.gz```
```bwa index -a bwtsw hg38.fa```
8. Once the index creation has finished, all the BWA index files are copied to */data/index* folder on all data nodes
9. Download and uncompress the evaluation datasets
```sh
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/
sequence_read/SRR062634_1.filt.fastq.gz
gunzip SRR062634_1.filt.fastq.gz

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12750/
sequence_read/ERR000589_1.filt.fastq.gz
gunzip ERR000589_1.filt.fastq.gz
```
10. Preprocess the FASTQ files into the custom MapReduce format by running the provided script
```fq_to_mrfastq.py <input_file.fastq>```
The output of this command is a text file called *output.mr.fastq*
11. Upload the custom preprocessed FASTQ files to HDFS by running the commands
```sh
hdfs dfs -mkdir -p /user/karl/mrbwa/output
hdfs dfs -put output.mr.fastq /user/karl/mrbwa/
```
12. Run MapReduce job using the command
```sh
yarn jar $HADOOP_HOME/hadoop-streaming-2.6.5.jar
-files app_mrbwa
-mapper "app_mrbwa/mapper_pybwa.py"
-combiner "app_mrbwa/combiner_pybwa.py"
-reducer "app_mrbwa/reducer_pybwa.py"
-input hdfs:///user/karl/mrbwa/output.mr.fastq
-output hdfs:///user/karl/mrbwa/output/alignment
```
13. Copy tsv output from HDFS to master node using command
```hdfs dfs -copyToLocal /user/karl/mrbwa/output/alignment/part-00000 .```

## Setting up MR-BWT-FM on Hadoop
It is assumed that an Apache Hadoop Cluster has already been setup and steps 1, 2, 3 from **MR-BWA** setup have already been done.
1. Upload all the provided Python scripts in *mr-bwt-fm* folder to the master node
2. Copy the uploaded Python scripts to *app_mrbwtfm* folder on the master node.
3. Assuming that human reference genome has already been downloaded in */data/index* folder on master node, the custom FM-Index is created by running 
```build_index.py <input_file> hg38_idx```
4. Once the index creation has finished, the hg38 idx index file is copied to */data/index* folder on all data nodes
5. Preprocess the FASTQ files downloaded in step 9 of  **MR-BWA** setup, into the custom MapReduce format by running the provided script ```parse_fq_file.py <input_file.fastq>```
The output of this command is a text file called *output.fq.reads*
6. Upload the custom preprocessed FASTQ files to HDFS by running the command
```sh
hdfs dfs -mkdir -p /user/karl/mrbwtfm/output
hdfs dfs -put output.fq.reads /user/karl/mrbwtfm/
```
7. Run MapReduce job using the command
```sh
yarn jar $HADOOP_HOME/hadoop-streaming-2.6.5.jar
-files app_mrbwtfm
-mapper "app_mrbwtfm/mapper_native.py"
-combiner "app_mrbwtfm/combiner_native.py"
-reducer "app_mrbwtfm/reducer_native.py"
-input hdfs:///user/karl/mrbwtfm/output.mr.reads
-output hdfs:///user/karl/mrbwtfm/output/alignment
```
8. Copy tsv output from HDFS to master node using command
```hdfs dfs -copyToLocal /user/karl/mrbwtfm/output/alignment/part-00000 .```

### Output Analysis
This provides an overview of the *output analysis.R* that is used to analyze the MapReduce output. The machine used to run this script should have [R](https://cran.r-project.org/doc/manuals/R-admin.html) installed. Running this script is straightforward:
1. Copy MapReduce output tsv file to some location on the laptop. *output analysis.R* assumes that output tsv file is called *mrbwa output.tsv*, however this can be easily changed.
2. Open laptop terminal and type R. R terminal should load up.
3. Run the command
```source(’output_analysis.R’)```
4. This will take some time to finish. When it’s done, it will create four bar charts based on the input tsv file:
-- *alignment_dist*: This shows the alignment operations distribution. The result shows base pairs count of matches, variations, insertions, deletions.
-- *variations_per_chromosome*: This shows the variations base pairs count grouped by chromosome.
-- *insertions_per_chromosome*: This shows the insertions base pairs count grouped by chromosome.
-- *deletions_per_chromosome*: This shows the deletions base pairs count grouped by chromosome.
