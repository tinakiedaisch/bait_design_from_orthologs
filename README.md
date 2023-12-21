# bait_design
Select target loci using orthologs


## 01. Input Data

This workflow assums you already have pruned orthologs from your taxa of interest.

- **Orthologs .tre files**

- **Homologs .tre file** (isoforms maseked)

- **Reference genome**
   - Annotations: 'your_reference.gff3' & Genomic sequence: 'your_reference.fa'

---
 ## 02. Extract gene names of your referece from the orthologs

> 💡 We need the the gene names (of the reference) in our orthologs so that we can later align the intron masked reference to the orthologs to split the exons.


This bash script gives you all gene number of your referenc in the output file 'output_gene_numbers.txt'

```bash
 #!/bin/bash

input_directory="/home/tree_files"
output_file="/home/tree_files/output_gene_numbers.txt"

for file in "$input_directory"/*.tre; do
    echo "Processing $file"
    result=$(grep 'your_reference' "$file" | sed -n 's/.*@\(AH[0-9]*\):.*/\1/p')
    if [ -n "$result" ]; then
        echo "$result" >> "$output_file"
    fi
done
```

---
## 03. Mask the Introns in Your Reference Genome

Follow the steps outlined in the [Evernote guide](https://www.evernote.com/shard/s383/client/snv?noteGuid=0a5af77a-8c6a-4ade-aab5-9a1044ac2817&noteKey=7d78636bd6b59130&sn=https%3A%2F%2Fwww.evernote.com%2Fshard%2Fs383%2Fsh%2F0a5af77a-8c6a-4ade-aab5-9a1044ac2817%2F7d78636bd6b59130&title=Prepping%2BZebrafinch%2BGenome%2Bfor%2BBlackbird%2Btranscript%2Balignments)


```bash
#create conda environment and installing BEDTools
micromamba create -n bedtools bedtools
```
```bash
#extract gene and Introns from the .gff3
 grep "gene" A_hyperchondriacus.gff3 >A_hyperchondriacus_genes.gff3
grep "exon" A_hyperchondriacus.gff3 >A_hyperchondriacus_exons.gff3
```
```bash
#subatract
bedtools subtract -a A_hyperchondriacus_genes.gff3 -b A_hyperchondriacus_exons.gff3 > A_hypochondriacus_introns.gff3
sed 's/gene/intron/g' A_hypochondriacus_introns.gff3 >A_hypochondriacus_introns_forreal.gff3
```
```bash
#Mask the intron sequences in the genome sequence using bedtools:
bedtools maskfasta -fi Ahypochondriacus_459_v2.0.fa -bed A_hypochondriacus_introns.gff3 -fo A_hypochondriacus_introns_masked.fa
```
```bash
#indexing
samtools faidx Ahypochondriacus_459_v2.0.fa
samtools faidx A_hypochondriacus_introns_masked.fa
```
```bash
#Extract gene sequences from the intron-masked genome.
bedtools getfasta -fi A_hypochondriacus_introns_masked.fa -bed A_hyperchondriacus_genes.gff3 -fo A_hypochondriacus_genes.fa -s
```
---
## 04. Pull Genes from Intron-Masked Genome

```bash
#!/bin/bash

awk '
  function is_header(line) {
    return (substr(line, 1, 1) == ">");
  }

  # Read the gene numbers from output_gene_numbers.txt into an array
  FNR == NR {
    gene_numbers[$1] = 1;
    next;
  }

  # Process the CDS file
  {
    # Check if the line is a header
    if (is_header($0)) {
      # Extract the gene number from the header
      gene_number = substr($0, 2);

      # Check if the gene number is in the list
      if (gene_number in gene_numbers) {
        # Print the header
        print $0;

        # Set a flag to print the sequence
        print_sequence = 1;
      } else {
        # Set a flag to skip the sequence
        print_sequence = 0;
      }
    } else if (print_sequence) {
      # Print the sequence
      print $0;
    }
  }
' output_gene_numbers.txt intron_masked_genome.fa > extracted_genes.fasta
```

💡 **Tip:** make sure that it is the same number of gene that you have extracted from the orthologs & check if they are the same

```bash
comm -23 sorted_list_extracted_genes.txt modified_sorted_output_gene_nr.txt > unique_to_sorted_extracted_genes.txt
### file should be empty
```

---
## 05. Create a seperate fasta file for each gene

```bash
awk '/^>/{if (filename) close(filename); filename=sprintf("gene_%s.fa", substr($0,2)); print > filename; next;} {print >> filename;}' extracted_genes.fasta
```

---
## 06. Write fastas from orthologous trees

```bash
cat *.fa > all.fa
```

To generate individual FASTA files for each tree using the Python script "write_fasta_from_trees," use the following command:

```bash
python2 ~/write_fasta_from_trees all.fa final_orthologs/
```

The Python script can be found [here](https://bitbucket.org/dfmoralesb/phylogenomic_dataset_construction/src/master/write_fasta_files_from_trees.py).

💡 make sure that all of them include your reference

Concatinate all FASTAS

```bash
cat *.fa > all.fa
```

---
## 07. Aligning fastas and add the intron masked reference to the alignmen

For a fast and simple multiple sequence alignment you can use MAFFT

```bash
#!/bin/bash

for file in *.fasta; do
    mafft --auto "$file" > "${file%.fasta}.aligned.fasta"
done
```
Then add you intron masked reference sequence to the alignmnet with the fuction mafft --add

```bash
#!/bin/bash

for file in *.fa; do
    mafft --add "$file" "${file%.fa}.aligned.fasta" > "${file%.fa}.add.fa"
done
```

---
## 08. Shorten sequence name to 8 characters

💡 Since we have to convert the fasta in phylip later, we have to unifiy and shorten the sequence name

```python
from Bio import SeqIO
import os

# Input directory
input_directory = "./"  # Adjust this to your input directory

# Iterate over all .add.fa files in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith(".add.fa"):
        input_file = os.path.join(input_directory, filename)
        output_file = os.path.join(input_directory, f"{filename[:-7]}_shortened.fa")

        # Open the input and output files
        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            # Iterate through the input FASTA file
            for record in SeqIO.parse(infile, "fasta"):
                # Shorten the sequence name to 8 characters
                shortened_name = record.id[:8]

                # Create a new SeqRecord with the shortened name
                modified_record = record
                modified_record.id = shortened_name
                modified_record.name = shortened_name
                modified_record.description = shortened_name

                # Write the modified record to the output file
                SeqIO.write(modified_record, outfile, "fasta")

```

---
## 09. Preparations for splitting Exons
### 09.1 converst FASTA to PHYLIP

💡we need the phylip format to extract colums in the next step

```bash
#!/bin/bash

for file in *_shortened.fa; do
    pxs2phy -s "$file" -o "${file%.add.fa}.phy"
done
```
### 09.2 mark introns with @

with this python script all clolums that have 'N' in the reference and gaps '-' in the other taxa will be replaced with '@'. This is done to make sure to split the complete exon in cases where the introns of the reference don't align with the other taxa.

```python
import sys
import glob

def process_file(input_filename):
    # Read the input file
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    # Extract the metadata line
    metadata = lines[0]

    # Extract the sequence lines, excluding the metadata line
    sequence_lines = lines[1:]

    # Get the length of the first sequence line
    first_sequence_length = len(sequence_lines[0].rstrip())

    # Check if all subsequent sequence lines have the same length as the first line
    if not all(len(line.rstrip()) == first_sequence_length for line in sequence_lines):
        print("Error: Sequence lines have different lengths.")
        return

    # Identify columns to replace
    columns_to_replace = set()
    for col in range(first_sequence_length):
        if sequence_lines[0][col] == '-' and (sequence_lines[-1][col].lower() == 'n' or sequence_lines[-1][col].lower() == 'n'):
            # Check if all intermediate rows have "-"
            intermediate_chars = [line[col] for line in sequence_lines[1:-1]]
            if all(char == '-' for char in intermediate_chars):
                columns_to_replace.add(col)

    # Replace identified columns
    output_lines = [metadata]  # Add back the metadata line to the output
    for line in sequence_lines:
        new_line = ''.join('@' if idx in columns_to_replace else char for idx, char in enumerate(line))
        output_lines.append(new_line)

    # Write the output to a new file
    output_filename = input_filename + ".replaced.phy"
    with open(output_filename, 'w') as file:
        file.writelines(output_lines)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide the directory path as an argument.")
    else:
        directory_path = sys.argv[1]

        # Find all files with the pattern _shortened.fa.phy in the specified directory
        file_pattern = f"{directory_path}/*_shortened.fa.phy"
        file_list = glob.glob(file_pattern)

        if not file_list:
            print(f"No files found with the pattern '{file_pattern}'.")
        else:
            for input_filename in file_list:
                process_file(input_filename)
                print(f"Processed: {input_filename}")
```

### 09.3 Convert PHYLIP back to FASTA

for downstream analysis we need again the FASTA format

```bash
#!/bin/bash

for file in *.fa.phy.replaced.phy; do
    pxs2fa -s "$file" -o "${file%.converted}.fa"
done
```

---
## 10. Splitting Exons

With this python script the exons wil be splitted according to the '@' signs in the alignments.

```python
def split_exons(fasta_file, outDIR):
    if not os.path.isabs(fasta_file):
        fasta_file = os.path.abspath(fasta_file)

    if os.path.isabs(outDIR) == False:
        outDIR = os.path.abspath(outDIR)
    if outDIR[-1] != "/":
        outDIR += "/"
    if outDIR == ".":
        outDIR = os.getcwd()

    # Create the output directory if it doesn't exist
    if not os.path.exists(outDIR):
        os.makedirs(outDIR)

    # Read the alignment
    align = AlignIO.read(fasta_file, "fasta")

    # Process the sequences
    for record in align:
        # Check if the sequence name starts with "AH000198"
        if record.name.startswith("AH000198"):
            exons = list(filter(None, re.split('@', str(record.seq))))
            number_exons = len(exons)
            introns = list(filter(None, re.split('a|c|t|g|-', str(record.seq))))
            number_introns = len(introns)
            number_modules = number_exons + number_introns
            coordinates = [0, len(exons[0])]

            if number_exons > 1:
                for i in range(0, number_exons-1):
                    coordinates += [coordinates[-1] + len(introns[i]), coordinates[-1] + len(introns[i]) + len(exons[i+1])]

            for i in range(0, number_exons):
                exon_align = align[:, coordinates[i*2]:coordinates[(i*2)+1]]
                file_name = outDIR + record.name + ".E" + str(i+1) + ".fna"

                # Write each exon to a separate file
                with open(file_name, "w") as output_handle:
                    AlignIO.write(exon_align, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python split_exons.py fasta_file outDIR")
        sys.exit(1)

    fasta_file, outDIR = sys.argv[1:]
    split_exons(fasta_file, outDIR)
```

---
## 11. Filtering exons

After the previous step you will get ten thousends of exons but which are probably too short, too uninformative, or contain gappy sequences. Therfor it is good to filter this one out first.

- Filter for length
- Filter for informativness (parsimony informative sites) with phykit [PhyKIT GitHub Repository](https://github.com/JLSteenwyk/PhyKIT)

```bash
#!/bin/bash

# Specify the input and output directories
input_dir="/path/to/your/input_directory"
output_dir="/path/to/your/output_directory"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all .fna files in the input directory
for file in "$input_dir"/*.fna; do
    # Extract the file name without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.fna}"

    # Run phykit command to calculate parsimony informative sites
    phykit pis "$file" > "$output_dir/${filename_no_ext}_pis.txt"

    echo "Parsimony Informative Sites for $filename calculated and saved to ${filename_no_ext}_pis.txt"
done
```

  💡 The output of the script will be displayed as follows:
col1: number of parsimony informative sites 
col2: total number of sites 
col3: percentage of parsimony informative sites

--> with this you can select a threshhold and remove exons below it


- Filter for missingness using trimal
  💡 Do not trimm the alignmnets. The following script will just list sequences with a proportion of gaps per sequence higher than 0.5. (50 % gaps). Then you can remove the one from the list.

```bash
   #!/bin/bash

# Set the path to your directory containing the .fna files
input_directory="/home/kiedaisch/kiedaisch_DATA_REMOTE/baits_design/final_orthologs/orthologs_AH/gene_files/Exons/uninformatives_removed/"

# Set the threshold value
threshold=0.5

# Iterate over all .fna files in the directory
for file in "$input_directory"*.fna; do
    # Extract the filename without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.fna}"

    # Run the script for each file
    /home/kiedaisch/apps/trimal/scripts/get_sequences_gaps_ratio.py -i "$file" --threshold "$threshold" > "$input_directory$filename_no_ext.trimal.txt"

    # Print a message indicating completion
    echo "Processed $filename"
done
```

remove sequences on the list with the following bash script:

```bash
for fasta_file in *.fna; do
    trim_file="${fasta_file%.fna}.trimal.txt"
    output_file="modified_${fasta_file}"

    # Check if the .trimal.txt file is not empty
    if [ -s "$trim_file" ]; then
        awk 'NR==FNR {sequencesToRemove[$2]; next} /^>/ {header=$0; sub(/^>/, "", header); if (!(header in sequencesToRemove)) print $0; next} {if (!(header in sequencesToRemove)) print}' "$trim_file" "$fasta_file" > "$output_file"
        echo "Processed $fasta_file"
    else
        # If .trimal.txt is empty, just copy the original .fna file
        cp "$fasta_file" "$output_file"
        echo "No modifications for $fasta_file"
    fi
done
```

---
## 11. Align Exons again with MACSE_OMM

<blockquote style="border-left: 4px solid #3498db; padding-left: 15px; margin-left: 0;">

Mafft did not align the sequences in frame, so if you are experiencing many frameshifts or badly aligned sequences, use `macse_omm`. This can be checked in Geneious Prime using the Translation option. The MACSE_OMM pipeline will align the sequences in respect to their AA translation, remove non-homologous sequence fragments, and trim extremities of the alignment.

</blockquote>


































 


