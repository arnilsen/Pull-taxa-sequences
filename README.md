# pull taxa sequences

pull_tax_sequences.py retrieves all sequences from a list of species and regions defined by the user from NCBI. The user supplies a txt file with one species per line and a list of regions/genes. The script downloads sequences for the regions and saves them into fasta files based on their regions. It also saves a csv file with all the species and regions available.


**System requirements**

pull_tax_sequences.py requires that python 3.x and biopython is installed. Installation instructions for python are available [here](https://www.python.org/downloads/), and biopython [here](https://biopython.org/wiki/Download).

Use git clone to obtain the script, or copy and paste the script into a text editor and save the file as pull_tax_sequences.py.


**Running the script**

```python
pull_tax_species.py -i species_list.txt -e your.email@random.com -r "internal transcribed spacer",rpb1,rpb2
```
Compulsory arguments are -i, -r, and -e. The list of regions/genes (-r) needs to be separated with commas. The email address is a requirement from NCBI so you can be contacted in case of excessive usage. The user can also specify the optional argument -v, --verbose. This will retain all temporary files. The script can take sometime to run, depending on the load at NCBI.



**Output**

pull_tax_species.py outputs two types of files, a csv and fasta files. The csv file lists species-vouchers with additional available markers. If no identifying voucher can be ascertained, then the species is listed with unknown appended to it. Several conveniences have been integrated into the csv file: identification of type specimens, identification of joint ITS and LSU sequences, and concatenation of gene names. The script can identify if a matching sequence is from a type specimen, appending ‘-TYPE’ to the marker. Additionally, if the sequence is from the ITS region and is longer than 700 bp long, the sequence is labelled as ‘ITS-LSU?’. Sequences that cannot be identified are binned as unknown. Sequences that contain multiple genes have their unique marker names concatenated e.g. trnL-trnF.
The fasta files contain all the sequences in the csv file with the same species-voucher identifier in the first position of the header. The second position contains the accession and third the marker name. This allows the user to use tools like grep to retrieve all makers from the generated fasta files.


**Additional**

If just the genus is supplied in the species list, the script will recover all species to a maximum of 5000 (set by retmax=5000 in the Entrez.esearch tool). The user is also warned if there are duplicate taxa or if their taxa are not available at NCBI. This script does not retrieve or parse whole genomes.


```
usage: pull_tax_species.py [-h] -i  -e  [-v] -r REGIONS [REGIONS ...]

Pull taxa sequences: pull all sequences for genes of interest from specified species

optional arguments:
  -h, --help            show this help message and exit
  -i , --input          enter file of species names, one species per line
  -e , --email          enter email address so NCBI can contact you if there
                        is a problem
  -v, --verbose         verbose, keep intermediate files
  -r REGIONS [REGIONS ...], --regions REGIONS [REGIONS ...]
                        regions to get, can be a comma seperated list. If
                        regions contain spaces it must be in "" marks

Additional information:
This script takes a file of species (one per line) and genes of interest and retrieves them from NCBI.
There are two main output files, a collection table and fasta files of the respective genes.
The collection table lists all the retrieved sequences with the accession numbers. The
collection names are formatted with the species name followed by their specimen voucher. This
naming convention is retained throughout to ensure seamless concatenation.
Several temporary files are also save to disk and deleted at the completion of the script. This
can be disabled with the --verbose flag.

pull_tax_species.py -i species_list.txt -e your.email@random.com -r "internal transcribed spacer",rpb1,rpb2

```
