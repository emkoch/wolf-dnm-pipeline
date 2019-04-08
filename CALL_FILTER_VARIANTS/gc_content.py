import pysam
import argparse
import os

def count_gc(ref_fasta, chrom, position, pad):
    sequence = ref_fasta.fetch(chrom, position-1-pad, position+pad)
    gc_count = 0
    for ii in range(len(sequence)):
        if sequence[ii] == "G" or sequence[ii] == "C":
            gc_count += 1
    return float(gc_count)

def remove_gc_range(vcf_file, outfile, stat_out, ref_fasta, check_perc, pad, max_gc_perc, min_gc_perc = 0.0):
    import os
    import sys
    import vcf
    import random
    import shutil
    import pysam
    import csv
    import gzip
    
    compressed = False
    if os.path.splitext(vcf_file)[1] == '.gz':
        compressed = True

    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]

    reference = pysam.FastaFile(ref_fasta)

    with open(outfile_write, 'w') as f_out:
        reader = vcf.Reader(filename=vcf_file, compressed=compressed)
        writer = vcf.Writer(stream = f_out, template = reader)
        record_count = 0
        kept_count = 0
        for record in reader:
            if random.random() <= check_perc :
                record_count += 1
                chrom = str(record.CHROM)
                position = int(record.POS)
                gc_prop = count_gc(reference, chrom, position, pad) / float(2 * pad + 1)
                if gc_prop <= max_gc_perc and gc_prop >= min_gc_perc:
                    kept_count += 1
                    writer.write_record(record)
        
    if os.path.splitext(outfile)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

    with open(stat_out, "w") as f:
        f.write(str(record_count) + "\n" + str(kept_count) + "\n")

    print("finished...")

    return record_count, kept_count
                
                

def gc_prop_vcf(vcf_file, ref_fasta, outfile, pad):
    import os
    import sys
    import vcf
    import random
    import shutil
    import pysam
    import csv
    import gzip
    
    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]

    compressed = False
    if os.path.splitext(vcf_file)[1] == '.gz':
        compressed = True
    
    reference = pysam.FastaFile(ref_fasta)
    
    with open(outfile_write, 'w') as f_out:
        outwriter = csv.writer(f_out, delimiter=',')
        outwriter.writerow(['chromosome', 'position', 'gc_prop'])
        
        reader = vcf.Reader(filename=vcf_file, compressed=compressed)
        for record in reader:
            chrom = str(record.CHROM)
            position = int(record.POS)
            gc_prop = count_gc(reference, chrom, position, pad) / float(2 * pad + 1)
            outwriter.writerow([chrom, position, gc_prop])

    if os.path.splitext(outfile)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_file", help="The VCF file to calculate GC content for")
    parser.add_argument("--reference", help="A fasta file for the reference genome")
    parser.add_argument("--outfile", help="The output csv to write")
    parser.add_argument("--pad", help="The number of sites around the focal site to include in calculation", type = int)
    args = parser.parse_args()
    
    gc_prop_vcf(args.vcf_file, args.reference, args.outfile, args.pad)
                        
    
if __name__ == "__main__":
    main()
    
