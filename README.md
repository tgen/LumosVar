# TumorOnly
Designed for calling somatic variants in tumor samples lacking matched normals.  Copy number and minor allele copy number are estimated by EM to determine the expected somatic and germline variant allele frequencies at each position.  Prior probabilities of germline or somatic variants at each position are calculated based on 1000 genomes allele frequencies or COSMIC counts.  A set of unmatched normals is used to find mean read depths and position quality scores.

PREREQUISITES
1. samtools and htslib (tested with 1.2)
http://www.htslib.org/download/
2. Matlab Runtime (MCR 9.0)
http://www.mathworks.com/products/compiler/mcr/

USAGE
Inputs are given in a config file in yaml format
