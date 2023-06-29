<br>

### Dealing with palindromic SNPs

<br>

When the reference strands of the exposure and outcome data sets are unknown, palindromic SNPs can be challenging to harmonise because the orientation of the effect allele is ambiguous. In some cases we can use the effect allele frequency to resolve the ambiguity. For example, consider a SNP with alleles A & T, with a frequency of 0.11 for allele A in the exposure data set and 0.91 in the outcome data set (see Table below). In both data sets the effect allele is A and assume they are from the same population and so are likely to have similar minor allele frequencies. Because the minor allele frequency in the exposure data set (0.11) corresponds to the A allele but the minor allele frequency in the outcome data (1-0.91 = 0.09) corresponds to the T allele, it is reasonable to infer that the exposure and outcome data are on different strands. Therefore, we need to reverse (multiply by -1) the outcome estimate to harmonise the estimates across data sets.

However, resolving these ambiguities is much more difficult when the effect allele frequencies are closer to 0.5. In such cases SNPs are determined to be too ambiguous to resolve, and are usually discarded. Note that this only applies when the exposure and outcome strands are uncertain. For more information see https://academic.oup.com/ije/article/45/6/1717/3072174?login=true

<br>

<p align="center">

<img src="https://i.ibb.co/SNgyRZQ/upload3.png" alt="Description" style="width:60%;"/>

</p>
<br>


The function has 3 options for dealing with palindromic SNPs, which are
identical to those used in the TwoSampleMR package:

(1) = assume SNPs in the exposure and outcome GWAS are on the same
    strand
(2) = try to infer the effect allele based on the effect allele
    frequency\*
(3) = Correct the strand for non-palindromic SNPs, but drop all
    palindromic SNPs

\*SNPs will be discarded if the exposure and outcome effect allele
frequency is \>0.42
