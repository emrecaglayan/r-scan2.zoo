import pysam
import sys

def extract_reads(bam_file, output_file, chrom, positions):
    # Open the BAM file
    bamfile = pysam.AlignmentFile(bam_file, "rb")

    # Handle both single and multiple positions
    if "," in positions:
        positions = [int(pos) for pos in positions.split(",")]
    else:
        positions = [int(positions)]

    # Dictionary to store read ID and the concatenated bases for each read
    read_bases = {}

    # Iterate through each position
    for pos in positions:
        # Fetch reads at the given position
        for pileupcolumn in bamfile.pileup(chrom, pos - 1, pos):  # zero-based position
            if pileupcolumn.pos == pos - 1:  # one-based to zero-based conversion
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # Get read ID and the base call at the current position
                        read_id = pileupread.alignment.query_name
                        base_call = pileupread.alignment.query_sequence[pileupread.query_position]

                        # Initialize or update the concatenated base call for each read
                        if read_id not in read_bases:
                            read_bases[read_id] = base_call
                        else:
                            read_bases[read_id] += base_call

    # Open a file to write the output
    with open(output_file, "w") as outfile:
        # Write the results to the output file
        for idx, (read_id, bases) in enumerate(read_bases.items(), start=1):
            outfile.write(f"{idx}\t{read_id}\t{bases}\n")

    # Close the BAM file
    bamfile.close()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_reads.py <bam_file> <output_file> <chromosome> <positions>")
        print("Example: python extract_reads.py input.bam output.txt chr1 100,200,300")
        sys.exit(1)

    bam_file = sys.argv[1]
    output_file = sys.argv[2]
    chrom = sys.argv[3]
    positions = sys.argv[4]

    extract_reads(bam_file, output_file, chrom, positions)

