![Scripting](https://img.shields.io/badge/Language-R-yellow.svg) ![Copyright](https://img.shields.io/badge/Copyright-(c)_2022_Max\_Stammnitz\_@TCG\_Cambridge-green.svg)

SubstitutionSafari – from substitution calls to substitution spectra
====================================================================
![SubstitutionSafari](/Images/SubstitutionSafari_logo.png)

Using the functions in [Substitution Safari](/SubstitutionSafari.R), we generate substitution – also termed single-base substitution (SBS) or single-nucleotide variant (SNV) – spectra in line with the widely-used classification scheme first adapted by [Alexandrov et al., 2013](https://www.nature.com/articles/nature12477) (see [COSMIC SBS signature catalogues](https://cancer.sanger.ac.uk/signatures/sbs/)).

Note that, in order to produce your own spectrum, you will need to specify:
* the reference genome fasta file based on which your alignments' substitution calls were generated
* a VCF or dataframe object containing one substitution per line, featuring the minimal set of columns "CHROM", "POS", "REF", "ALT"
* a table with the triplet frequency counts of your reference genome or exome, which can be used to normalise the observed substitution counts for context-specific mutation opportunity (an example table can be found [here](/Files/reference_genome_trinucleotides.txt))

My code then groups all of your substitutions with reference to the 96 different, pyrimidine-centered subgroups of C>A, C>G, C>T, T>A, T>C, and T>G types. Plotted spectra look like this example:

![example](/Images/Example_spectrum.png)

In the above case, we see enrichments of C>T substitutions, particularly in A[C]G and G[C]G contexts. These spikes indicate the presence of COSMIC signature [SBS6](https://cancer.sanger.ac.uk/signatures/sbs/sbs6/), a mutational signature commonly identified in tumours with DNA mismatch repair deficiency. You might want to corroborate this by inspecting and fitting combinations of COSMIC signatures to the spectrum, for example by using a tool such as [sigfit](https://github.com/kgori/sigfit).

Note that the script and its associated results are yet unpublished/peer-reviewed. It has nevertheless already been extensively benchmarked with substitution calls made against the Tasmanian devil reference. In theory, this code should run smoothly for ANY species with a reference genome. 

Find it helpful or you simply enjoy making substitution spectra against your brand-new species reference genome with these lines? Then why not follow the same style with your indel calls – have a look at its sister library [Indelwald](https://github.com/MaximilianStammnitz/Indelwald)! If you do face a challenge in using the code or wish to provide general feedback, please get in touch directly via maxrupsta@gmail.com
