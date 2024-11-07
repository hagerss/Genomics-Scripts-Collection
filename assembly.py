"""This program is designed to perform de novo assembly of sequencing reads directly
 from a file, employing a greedy algorithmic approach. It requires a minimum overlap of
 k=10 nucleotides to merge reads into contigs, ensuring a non-divergent assembly process.
After assembling the reads into contigs, the program maps the original reads back to the contigs
 to calculate and visualize coverage along each contig."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import matplotlib.pyplot as plt

def findOverlap(a, b, k):
    """Find initial overlap position between two sequences.

    :param
    a : The first sequence.
    b (str): The second sequence to compare against the first.
    k (int): Minimum number of nucleotides for a valid overlap.

    Returns:
     The starting index of the overlap in sequence a, -2 for invalid overlap, or -1 for no overlap.
    """
    findSeq = b[0:k]  # Get the initial sequence of 'b' of length 'k'
    lp = a.find(findSeq)  # Look for this sequence in 'a'
    if lp + len(b) > len(a):  # Check if 'b' extends beyond 'a' indicating a valid overlap
        return lp
    elif lp > 0:  # Found overlap, but it's not valid for merging
        return -2
    else:  # No overlap found
        return -1


def mergeOverlaps(reads, k):
    """Merge reads with an overlap of at least 'k' nucleotides.

    :param:
    reads : The list of sequence reads.
    k  Minimum overlap length.

    Returns:
    list: List of merged sequence reads.
    """
    # Initial prints for debugging
    print("test")  # Debugging print
    print(reads[0])  # Show the first read for context
    for i in range(len(reads)):  # Iterate through reads
        for j in range(len(reads)):  # Compare each read with every other
            if (len(reads[j]) > k) and (len(reads[i]) > k) and (
                    i != j):  # Ensure both reads are long enough and not the same
                ovlp = findOverlap(reads[i], reads[j], k)  # Check for overlap
            else:
                ovlp = -1  # Set overlap to -1 if conditions are not met
            if ovlp > -1:  # If valid overlap found
                reads[i] = reads[i][0:ovlp] + reads[j]  # Merge reads
                reads[j] = ""  # Empty merged read
            if ovlp == -2:  # If overlap found but not valid
                reads[j] = ""  # Remove read
    return [read for read in reads if read]  # Filter out empty strings


def getSeq(filename):
    """Read sequences from a file.

    :param:
    filename (str): Path to the file containing sequences.

    Returns:
    list: List of sequence reads.
    """
    with open(filename, 'r') as f:  # Open file
        allreads = [line.rstrip() for line in f]  # Read and strip newline characters
    return allreads


def assembleReads(filename, k):
    """Assemble reads into contigs with a minimum overlap of 'k' nucleotides.

    Args:
    filename (str): Path to the file containing sequences.
    k (int): Minimum overlap length.

    Returns:
    list: List of assembled contigs.
    """
    testread = getSeq(filename)  # Get sequences from file
    xP = len(testread) + 1  # Set previous length
    x = len(testread)  # Set current length
    # While the length of testread decreases
    while (xP > x):
        xP = x  # Update previous length
        testread = list(filter(None, mergeOverlaps(testread, k)))  # Merge overlaps and remove empty reads
        x = len(testread)  # Update current length
        print(x, "   ", xP)  # Debugging print
    cNum = 0  # Initialize contig counter

    return testread  # Return assembled contigs


def writeContigsToFasta(contigs, filename="assembled_contigs.fasta"):
    """Writes assembled contigs to a FASTA file using BioPython.

    Args:
    contigs (list): List of assembled contigs.
    filename (str): Output FASTA file name.
    """
    seq_records = []
    for i, contig in enumerate(contigs, 1):
        seq_record = SeqRecord(Seq(contig), id=f"contig_{i}", description="")
        seq_records.append(seq_record)
    SeqIO.write(seq_records, filename, "fasta")


def calculateCoverage(contigs, reads):
    """
    Calculates coverage for each contig based on the list of reads.

    :param:
    contigs : dictionary of the assembled contigs.
    reads : The original sequence reads.

    Returns:
    dict: A dictionary with contig names as keys and lists of coverage values as values.
    """
    coverage = {}  # Initialize an empty dictionary for coverage data
    for contig_name, contig_sequence in contigs.items():  # Iterate through each contig
        coverage[contig_name] = [0] * len(contig_sequence)  # Initialize coverage list with zeros for each base in contig
        for read in reads:  # Iterate through each read
            start = contig_sequence.find(read)  # Find the start position of the read in the contig
            while start != -1:  # If the read is found within the contig
                for i in range(len(read)):  # Iterate through each base in the read
                    coverage[contig_name][start + i] += 1  # Increment coverage for each base of the read in the contig
                start = contig_sequence.find(read, start + 1)  # Look for next occurrence of the read in the contig
    return coverage  # Return the coverage data


def plotCoverage(coverage, output_file="coverage_plot.png"):
    """Plots coverage for each contig and saves the plot to a file.

    Args:
    coverage (dict): Coverage data.
    output_file : Filename for the output plot.
    """
    plt.figure(figsize=(10, 6))
    for contig_name, contig_coverage in coverage.items():
        plt.plot(contig_coverage, label=contig_name)
    plt.xlabel("Position in Contig")
    plt.ylabel("Coverage")
    plt.title("Coverage Along Contigs")
    plt.legend()
    plt.savefig(output_file)
    plt.close()


def main():

    filename = "seqReadFile20"
    k = 10
    assembled_contigs = assembleReads(filename, k)
    contig_dict = {f"contig_{i + 1}": contig for i, contig in enumerate(assembled_contigs)}
    reads = getSeq(filename)
    # Calculate and plot coverage
    coverage = calculateCoverage(contig_dict, reads)
    plotCoverage(coverage, "coverage_plot.png")
    writeContigsToFasta(assembled_contigs, "assembled_contigs.fasta")
    print("Assembly and coverage analysis complete. Outputs have been saved.")

if __name__ == "__main__":
    main()

