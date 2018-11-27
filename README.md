# Mitochondrial-variation-detection
Scripts used in article 'A simple pipeline for heteroplasmic variation detection of the mitochondrial genome from whole genome sequencing data'

Use `$ git clone git@github.com:duanmq1994/Mitochondrial-variation-detection.git` download this source codes or directly download the ZIP.
     

### For the heteroplasmic variation detection
This step suggests that the users align the sequencing data against two references(**MT.fasta** and **MT-8k.fatsa**) respectively.     
Then the .bam/.sam files need to be kept and used in the following steps.     
Or users can creat their own .fasta files using the **fa8k.pl**.      
This script will creat a new .fasta file of mitochondrial genome to eliminate the influence on the head and tail reads because of a circular reference.      

- Use the base_indel_count.pl.      
This script will read the .bam/.sam file (use normal reference and 8k reference respectively) which has been sorted by SAMtools, then output three file: **base-count.tsv**, **indel-count.tsv** and **point-variation.tsv(alternative)**.

- Use the **1_8k_merge.pl** to merge those two sets of .tsv files.


### For the linkage detection     
This step also needs users to creat the following results files using .bam/.sam files from two reference.

- Use the **indel_base_link.pl**.
This script will use the .bam/.sam file which aligned against two reference from *heteroplasmic variation detection* step. It will creat a primary **id-variation.tsv** file (intermediate result).     
- Use the **variations_merge.pl** to get the **linkages-count.tsv** file.      
- Use the **linkage_classify.pl**.    
This script will read the **linkages-count.tsv** file and output three kinds of linkage files: **indels-linkage.tsv**, **indel-point_variations-linkage.tsv** and **point-variations.tsv**.
