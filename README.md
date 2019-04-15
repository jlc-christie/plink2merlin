# plink2merlin
A simple command line utility written in Python to convert from plink format to Abecasis Lab's QTDT format (used by Merlin). Although similar, plink's .ped & .map are not directly compatible with the QTDT format. This also produces a compatible .map file to be used by specific functions of Merlin (e.g. IBD calculations). 

This script still needs work, feel free to submit pull requests. There is a function call that has been commented out in the script which filters families with >1 offspring, this was used to extract only families with valid sib-pairs, if this is useful and are comfortable editing scripts, just uncomment the line in the `process_chrom` function.

## Usage:

```
plink2merlin.py -h
usage: plink2merlin.py [-h] [--bed BED] [--ped PED] [--vcf VCF]
                       [--plink-binary PLINK_BINARY]

File format converter to take a plink binary file or vcf and convert it to a
format which MERLIN accepts, which consists of a .ped and a .dat file. It is
important to note that the .ped file here is not the same as the plink .ped
file.

optional arguments:
  -h, --help            show this help message and exit
  --bed BED             PLINK binary input file path stub
  --ped PED             PLINK input file path stub
  --vcf VCF             VCF file path
  --plink-binary PLINK_BINARY
                        Path to PLINK binary, useful if plink isn't globally
                        accessible or you want to use a specific version of
                        plink
```

## Examples
The only required argument of the utility is the input file, of which the format must be specified. Despite the usage information, the only currently supported file format is plink's binary (.bed,.bim,.fam) format. Although, a lot of the work for supporting plink's text format (.ped, .map) and .vcf is already done (feel free to contribute to finish it off). 

### Basic usage
Below is an example of running the script, assuming `plink_binary_stub.bed`, `plink_binary_stub.bim` & `plink_binary_stub.fam` are all in the same directory as the script.
```bash
python plink2merlin.py --bed plink_binary_stub
```

### Supplying a plink binary
The utility allows for passing the plink binary directly, in case plink is installed locally, or under a different alias, here the plink binary is in the same directory as the script.
```bash
python plink2merlin.py --plink-binary ./plink --bed genome_data
```

Another example where the binary is passed as an absolute path:
```bash
python plink2merlin.py --plink-binary /home/my-username/plink_binaries/plink1.9 --bed genome_data
```

## Output 
Output will be saved in to a created `merlin_input_files` directory, named by chromosome, this can all be changed by editing the script but not by adding command line flags (again, feel free to submit a PR if you implement this). 

## How it works 
1. Split data by chromosome using plink, save to temporary directory (merlin does not account for chromosomes)
2. For each chromosome:
    1. Copy data from .fam in to data structures in memory for fast and easy access.
    2. Fill in missing parent data (plink accepts parent IDs which aren't in the dataset, merlin doesn't, so we have to add dummy entries with null genotypes)
    3. Create dict mapping family ID to list of the individuals in the family, keep in memory
    4. Find and split disjoint families (plink accepts disjoint families (families under the same FID that aren't _explicitly_ related), merlin throws exception if it finds them)
    5. Get genotype/variant information
        1. Prune SNP data using plink so that 1.) SNPs aren't in LD, and 2.) uncommon variants aren't included 
        2. Create plink `.ped`/`.map` files from binary, including only pruned SNPs
        3. extract genotype data from `.ped` and variant data from `.map`
        4. Create dict mapping from fid+iid as key and genotype data as value
    6. Using all information gathered, parse data in to QTDT records, create directory for ouput and save `.ped`, `.map` & `.dat` to file
3. Delete temporary files and cleanup
