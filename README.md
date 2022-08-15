# VSensus
VSensus is a simple consensus caller that was designed for viral samples due to their ultra-deep nature and higher mutaiton rate. When running VSensus you can vary the coverage threshold with which to call a consensus (sites below this will yield an N), vary the quality threshold for individual bases to be utilised in the consensus calling, and vary the variant frequency with which to use an ambiguity code.

Vsensus works off the [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) format generated from [samtools](https://samtools.github.io).

To create an mpileup file from a BAM file you need to use the 'samtools mpileup' command. As viral samples are ultra deep we recommned over-riding some of the mpileup defaults to produce an mpileup file with all the data:

```
samtools mpileup -aa -d 0 -A -Q 0 -f ref.fasta sample.bam > sample_mpileup.txt
```

* -aa = output absolutely all positions, including unused (i.e. zero coverage) reference positions
* -d 0 = max per-file depth - the default is 8000 - setting to 0 removes any threshold on depth (we want to use all the data to call the consensus)
* -A = do not discard anomalous read pairs
* -Q 0 = skip bases with quality less than 0 (so ouput all bases). The mpileup default is 13. VSensus can filter bases below a threshold downstream, so we recommned oututting all the data into mpileup file so it only needs to be created once.
* -f ref.fasta = the reference sequence file used to align the reads to when the BAM file was created
* sample.bam = the input aligbment BAM file to process
* sample_mpileup.txt = direct the mpileup output into a text file

Once an mpileup file has been created VSensus can be run to create an consensus FASTA sequence. VSensus is written in [Java](https://www.java.com/en/) which means the jar file file can be executed on any machine with a java runtime environment (JRE) version 8 (1.8) and above (although it has currently only been tested on Mac OS-X and Ubuntu Linux).

To run VSensus simply type this into the command line (where FOLDER/PATH is the path to the location of the VSENSUS.jar file)

```
java -jar /FOLDER/PATH/VSENSUS.jar sample_mpileup.txt
```

This will use the defaults of coverage threshold = 5 (-c=5), minimum qulaity = 0 (-q=0) and create an output file called INPUT_vsensus.fasta (so in the above it would be sample_mpileup_vsensus.fasta). While it is running, it will output messages to the terminal about mutations accepted into the consensus, and warnings about high frequency mutations and indels that didn't quite make it to the consensus level.

The defaults can be changed via the command line by providing arguments after the mpileup filename:

```
java -jar /FOLDER/PATH/VSENSUS.jar sample_mpileup.txt -c=10 -q=20
```

To vary the threshold at which an ambiguity code is used, use the -a option:
* -a=ambiguityFreq : [double] the secondary nucleotide frequency needed to trigger Ambiguity codes in the consensus: default >= 0.4 (i.e. the 2nd highest base at the position must have a frequency >40% for a ambibuity code to be used, an alternate view/explanation to this is that the top base will typically need to have a frequency of >60% to avoid an ambiguity code being used).

There are a number of other options currently in testing, a manuscript is under preparation along with a large number of abnormal/divergent virus samples which trigger unexpected behaviour in other variant callers.


