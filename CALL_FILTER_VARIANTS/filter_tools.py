import itertools
import shlex
import subprocess
import sys
import gzip
import os
import csv
import codecs
from platform import python_version

def get_homo_pls(infiles, outfile, child, mother, father):
    import vcf
    import os
    import csv
    
    with open(outfile, "w") as out_f:
        out_writer = csv.writer(out_f, delimiter=",")

        out_writer.writerow(["CHR", "POS", "REF_MOTHER", "REF_FATHER", 
                             "PL_0_0_MOTHER", "PL_0_1_MOTHER", "PL_1_1_MOTHER", 
                             "PL_0_0_FATHER", "PL_0_1_FATHER", "PL_1_1_FATHER"])
        
        for infile in infiles:
            reader = vcf.Reader(filename = infile, 
                                compressed = os.path.splitext(infile)[1] == ".gz")
            for record in reader:
                if(record.genotype(child)["GT"] == "0/0" and 
                   record.genotype(mother)["GT"] == "0/0" and 
                   record.genotype(father)["GT"] == "0/0"):
                    out_writer.writerow([record.CHROM, record.POS,
                                         record.genotype(mother)["AD"][0], record.genotype(father)["AD"][0]] + 
                                        record.genotype(mother)["PL"] + record.genotype(father)["PL"])

def get_trans_pls(infiles, outfile, child, mother, father, both = False):
    import vcf
    import os
    import csv
    
    with open(outfile, "w") as out_f:
        out_writer = csv.writer(out_f, delimiter=",")
        
        out_writer.writerow(["CHR", "POS", "DP_CHILD", "REF_CHILD", "ALT_CHILD",
                             "DP_MOTHER", "REF_MOTHER", "ALT_MOTHER",
                             "DP_FATHER", "REF_FATHER", "ALT_FATHER",
                             "PL_0_0", "PL_0_1", "PL_1_1"])

        for infile in infiles:
            reader = vcf.Reader(filename = infile, 
                                compressed = os.path.splitext(infile)[1] == ".gz")
            for record in reader:
                if both:
                    if((record.genotype(mother)["GT"] == "0/1" and 
                        record.genotype(father)["GT"] == "0/1") and 
                       record.genotype(child)["GT"] == "0/1"):
                        out_writer.writerow([record.CHROM, record.POS, 
                                             record.genotype(child).data.DP, 
                                             record.genotype(child)["AD"][0], record.genotype(child)["AD"][1],
                                             record.genotype(mother).data.DP, 
                                             record.genotype(mother)["AD"][0], record.genotype(mother)["AD"][1],
                                             record.genotype(father).data.DP, 
                                             record.genotype(father)["AD"][0], record.genotype(father)["AD"][1],
                                             record.genotype(child)["PL"][0], record.genotype(child)["PL"][1],
                                             record.genotype(child)["PL"][2]])
                else:
                    if((record.genotype(mother)["GT"] == "0/1" or 
                        record.genotype(father)["GT"] == "0/1") and 
                       record.genotype(child)["GT"] == "0/1"):
                        out_writer.writerow([record.CHROM, record.POS, record.genotype(child).data.DP,
                                             record.genotype(child)["AD"][0], record.genotype(child)["AD"][1],
                                             record.genotype(mother).data.DP, 
                                             record.genotype(mother)["AD"][0], record.genotype(mother)["AD"][1],
                                             record.genotype(father).data.DP, 
                                             record.genotype(father)["AD"][0], record.genotype(father)["AD"][1],
                                             record.genotype(child)["PL"][0], record.genotype(child)["PL"][1],
                                             record.genotype(child)["PL"][2]])
                    
    

def deslowmo_vcf_to_csv(infiles, outfile, child, mother, father, stats):
    import vcf
    import shutil
    import csv
    import os

    with open(outfile, "w") as out_f:
        out_writer = csv.writer(out_f, delimiter=",")
        
        out_writer.writerow(["CHR", "POS", "ALT_CHILD", "DP_CHILD", "SUM_CHILD",
                             "DP_MOTHER","DP_FATHER", "PP", "call_child", "QUAL"] + stats)
        
        for infile in infiles:
            reader = vcf.Reader(filename = infile, 
                                compressed = os.path.splitext(infile)[1] == ".gz")
            for record in reader:
                out_writer.writerow([record.CHROM, record.POS, record.genotype(child).data.AD[1],
                                     record.genotype(child).data.DP, sum(record.genotype(child).data.AD),
                                     record.genotype(mother).data.AD[0],
                                     record.genotype(father).data.AD[0], 
                                     record.INFO["pp"],
                                     record.genotype(child)["GT"], record.QUAL] + 
                                    [record.INFO[stat] for stat in stats])
                

def pp_cutoff(infile, outfile, cutoff):
    import sys
    import csv
    import gzip
    import collections
    import os
    import vcf
    import subprocess
    import math
    import itertools
    import numpy as np
    cutoff = float(cutoff)
    
    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]
        
    with open(outfile_write, "w") as fout:
        reader = vcf.Reader(filename = infile, 
                            compressed = os.path.splitext(infile)[1] == ".gz")
        writer = vcf.Writer(stream=fout, template=reader)
        for record in reader:
            if record.INFO["pp"] > cutoff:
                writer.write_record(record)

    if os.path.splitext(outfile)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

def deslowmo(infile, outfile_vcf, outfile_pp, child, mother, father, prior, theta, stats):
    import vcf
    import shutil
    import csv
    import os
    import collections
    import itertools
    import codecs
    import numpy as np
    from deslowmo_utils import index, trans_matrix, get_sample, check_record

    theta = float(theta)

    t_mat = trans_matrix(prior)
    
    outfile_write = outfile_vcf
    if os.path.splitext(outfile_vcf)[1] == '.gz':
        outfile_write = os.path.splitext(outfile_vcf)[0]

    with open(outfile_write, "w") as f_out, open(outfile_pp, "w") as f_out_pp:
        reader = vcf.Reader(filename=infile, compressed=True)
        # Edit info fields before opening writer using this as a template
        Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
        pp_field = Info('pp', 1, 'Float', 'posterior probability', 'denovogear model', None)
        reader.infos['pp'] = pp_field

        writer = vcf.Writer(stream=f_out, template=reader)
        pp_writer = csv.writer(f_out_pp, delimiter=",")
        pp_writer.writerow(["CHR", "POS", "ALT_CHILD", 
                            "DP_" + mother,"DP_" + father, "PP"] + stats)
        
        for record in reader:
            pp_out = calc_pp(record, child, mother, father, t_mat, theta)
            pp_total = pp_out[1] + pp_out[25]
            record.INFO['pp'] = pp_total
            writer.write_record(record)
            stat_values = []
            for stat in stats:
                if stat in record.INFO.keys():
                    stat_values.append(record.INFO[stat])
                else:
                    stat_values.append("NA")
            pp_writer.writerow([record.CHROM, record.POS, record.genotype(child).data.AD[1],
                                record.genotype(mother).data.AD[0],
                                record.genotype(father).data.AD[0],
                                pp_total] + stat_values)
            
    if os.path.splitext(outfile_vcf)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

        
def calc_pp(record, c, m, f, tm, theta):
    import numpy as np
    from deslowmo_utils import index, trans_matrix, get_sample, check_record
    result = np.array([0.0]*27)
    lh_c = get_sample(record, c)["PL"]
    lh_m = get_sample(record, m)["PL"]
    lh_f = get_sample(record, f)["PL"]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ind = index(i,j,k)
                prhet = ( (j!=1)*(1-theta) + (j==1)*theta )*( (k!=1)*(1-theta) + (k==1)*theta )
                result[ind] = -lh_c[i]/10.0 - lh_m[j]/10.0 - lh_f[k]/10.0 + tm[ind] + np.log10(prhet)
    ## Use log-sum-exp trick to calc probabilities
    b = np.max(result)
    logsum = b + np.log10( np.sum( np.power(10.0,result - b) ) )
    return np.power(10.0, (result - logsum))


def count_variant_bins(infile, outfile, binwidth=400):
    import vcf
    with open(infile, "r") as f_in, open(outfile, "w") as f_out:
        reader = vcf.Reader(f_in)
        # initialize the current bin 
        current_bin = [1, 100]
        current_count = 0
        for record in reader:
            pos = record.POS
            # if the current record is in the current bin add to count
            if current_bin[0] <= pos <= current_bin[1]:
                current_count += 1
            # if current record is not in the current bin, write the current bin count
            # and move the bin to the current record
            else:
                f_out.write(str(current_count) + "\n")
                current_count = 0
                while not (current_bin[0] <= pos <= current_bin[1]):
                    current_bin[0] += 100
                    current_bin[1] += 100
                    # if in new bin add to count
                    if current_bin[0] <= pos <= current_bin[1]:
                        current_count += 1
                        # if record not in new bin set count to zero and write this
                    else:
                        current_count = 0
                        f_out.write(str(current_count) + "\n")
        f_out.write(str(current_count) + "\n")

def remove_parental_nonref(infile, outfile, parents):
    import vcf
    '''
    Write all sites from infile to outfile where the parents were
    not CALLED as having any alt alleles
    '''
    with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
        reader = vcf.Reader(f_in)
        writer = vcf.Writer(stream = f_out, template = reader)
        total_count = 0
        kept_count = 0
        for record in reader:
            total_count += 1
            has_nonref = False
            for samp in record.samples:
                if samp.sample in parents:
                    if samp.gt_type != 0:
                        has_nonref = True
            if not has_nonref:
                writer.write_record(record)
                kept_count += 1
        writer.close()
    print('Kept {:d} of {:d} positions for a percentage of {:f}\%'.format(kept_count, 
                                                                          total_count, 
                                                                          kept_count/
                                                                          float(total_count)))
    return (total_count, kept_count)

def remove_high_depth(infile, outfile, depthcut):
    import vcf
    '''
    Write all sites from infile to outfile where no individuals in
    the vcf has a coverage greater than depthcut
    '''
    with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
        reader = vcf.Reader(f_in)
        writer = vcf.Writer(stream = f_out, template = reader)
        total_count = 0
        kept_count = 0
        for record in reader:
            total_count += 1
            depths = get_depths(record)
            if sum([depth >= depthcut for depth in depths]) == 0:
                writer.write_record(record)
                kept_count += 1
        writer.close()
    print('Kept {:d} of {:d} positions for a percentage of {:f}\%'.format(kept_count, 
                                                                          total_count, 
                                                                          kept_count/
                                                                          float(total_count)))
    return (total_count, kept_count)

def remove_depth_range(infile, outfile, min_depth, max_depth):
    '''
    Write all sites from infile to outfile where no individuals in
    the vcf has a coverage greater than depthcut
    '''
    import os
    import vcf
    import gzip
    import shutil

    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]

    with open(outfile_write, 'w') as f_out:
        reader = vcf.Reader(filename=infile, compressed=True)
        writer = vcf.Writer(stream = f_out, template = reader)
        total_count = 0
        kept_count = 0
        for record in reader:
            total_count += 1
            depths = get_depths(record)
            # Require that all depths are within the specified range for 
            # all individuals in the sample
            if (sum([depth > max_depth for depth in depths]) == 0 and 
                sum([depth < min_depth for depth in depths]) == 0):
                writer.write_record(record)
                kept_count += 1
        writer.close()

    if os.path.splitext(outfile)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

    print('Kept {:d} of {:d} positions for a percentage of {:f}\%'.format(kept_count, 
                                                                          total_count, 
                                                                          kept_count/
                                                                          float(total_count)))
    return (total_count, kept_count)


def get_potential_denovos(infile, outfile, child):
    import vcf
    import gzip
    import shutil
    
    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]
        
    with open(outfile_write, "w") as f_out:
        reader_in = vcf.Reader(filename=infile, compressed=True)
        writer = vcf.Writer(stream = f_out, template = reader_in)
        
        for record in reader_in:
            if 'AD' in record.genotype(child).data._fields:
                if record.genotype(child).data.AD[1] > 0:
                    writer.write_record(record)
                    
        writer.close()

    if os.path.splitext(outfile)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

    return None

def cal_coverage_histogram(infile, outfile, depthcut):
    import vcf
    import os
    import csv
    
    assert os.path.splitext(infile)[1] == ".gz"

    reader_in = vcf.Reader(filename=infile, compressed=True)
    
    wrote_header = False

    with open(outfile, "w") as csvfile:
        hist_writer = csv.writer(csvfile, delimiter=",")
        for record in reader_in:
            if not wrote_header:
                sample_names = [record.samples[ii].sample for ii in range(len(record.samples))]
                hist_writer.writerow(["name"] + sample_names)
                wrote_header = True
                sample_bins = dict()
                for sample in sample_names:
                    sample_bins[sample] = [0]*(depthcut + 1)
            for sample in sample_names:
                depth = record.genotype(sample)["DP"]
                if depth is None:
                    sample_bins[sample][0] += 1
                elif depth >= depthcut:
                    sample_bins[sample][depthcut] += 1
                else:
                    sample_bins[sample][depth] += 1
        for ii in range(depthcut + 1):
            hist_writer.writerow([ii] + [sample_bins[sample][ii] for sample in sample_names])
        
    return None

def remove_alt_read_parental_vcf(infile, outfile, mother, father, readcut):
    import os
    import vcf
    import gzip
    import shutil

    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]

    with open(outfile_write, "w") as f_out:
        reader_in = vcf.Reader(filename=infile, compressed=True)
        writer = vcf.Writer(stream = f_out, template = reader_in)
        total_count = 0
        kept_count = 0
        for record in reader_in:
            total_count += 1
            passes_count = True
            if 'AD' in record.genotype(mother).data._fields:
                if record.genotype(mother).data.AD[1] > readcut:
                    passes_count = False
            if 'AD' in record.genotype(father).data._fields:
                if record.genotype(father).data.AD[1] > readcut:
                    passes_count = False
            if passes_count:
                kept_count += 1
                writer.write_record(record)
        writer.close()

    if os.path.splitext(outfile)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

    print('Kept {:d} of {:d} positions for a percentage of {:f}\%'.format(kept_count, 
                                                                          total_count, 
                                                                          kept_count/
                                                                          float(total_count)))
    return (total_count, kept_count)

def iter_count_flank_variants(chrom, pos, varfile, flank):
    import vcf
    import os
    window = [pos - flank, pos + flank +1]
    count = 0
    with open(varfile, "r") as ff:
        reader = vcf.Reader(ff)
        for record in reader:
            if window[0] <= record.POS <= window[1] :
                count += 1
            elif record.POS > window[1]:
                break
    return count

def fast_remove_high_var_gzip(infile, varfile, outfile, varcut, window):
    print("new version!!")
    import vcf
    import copy
    '''
    **FAST VERSION**
    Write all sites from infile to outfile where no individuals in
    the vcf have a number of variants as great or greater than varcut
    in the varfile within the specified window size.
    - The window specifies the number of bases on either side of the variant position to consider
    - TO BE CLEAR, THE PUTATIVE DENOVO VARIANT IS INCLUDED IN THIS CALC
    '''
    total_count = 0
    kept_count = 0
    with open(outfile, "w") as f_out:
        reader_in = vcf.Reader(filename=infile, compressed=True)
        reader_var = vcf.Reader(filename=varfile, compressed=True)
        writer = vcf.Writer(stream=f_out, template=reader_in)
        # initialize position of variant reader and previous region variant list
        try:
            current_var = next(reader_var)
        except StopIteration:
            print("no variants in the file...")
            return (0, 0)
        prev_varlist = []
        for record in reader_in:
            total_count += 0
            current_varlist = []
            window_bottom = int(record.POS) - window
            window_top = int(record.POS) + window
            # check previous variants and add those in region to the current region variant list
            current_pos_is_var = False
            for var in prev_varlist:
                if window_bottom <= var.POS <= window_top :
                    if var.POS == record.POS:
                        current_pos_is_var = True
                    current_varlist.append( copy.copy(var) )
            while current_var.POS <= window_top :
                if current_var.POS == record.POS:
                    current_pos_is_var = True
                if window_bottom <= current_var.POS <= window_top :
                    current_varlist.append( copy.copy(current_var) )
                try:
                    current_var = next(reader_var)
                except StopIteration:
                    print("reached file's end...")
                    break
            var_count = len(current_varlist)
            if current_pos_is_var:
                var_count -= 1
            if var_count < varcut:
                writer.write_record(record)
                kept_count += 1
            prev_varlist = copy.deepcopy(current_varlist)
    return (total_count, kept_count)

def select_snps_only(infile, outfile):
    """
    Go through entries in the input vcf and only output snp sites
    """
    import os
    import vcf
    import gzip
    import shutil
    
    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]
    
    with open(outfile_write, 'w') as f_out:
        reader_in = vcf.Reader(filename=infile, compressed=True)
        writer = vcf.Writer(stream=f_out, template=reader_in)
        for record in reader_in:
            if record.is_snp:
                writer.write_record(record)
    
    if os.path.splitext(outfile)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

    return None

def get_stats(varfiles, outfile, stats):
    """
    Looks through the info field in all varfiles for each stat in stats and
    writes that to a csv outfile
    """
    import vcf
    import csv

    with open(outfile, "w") as f_out:
        outwriter = csv.writer(f_out, delimiter=",")
        # write the heading
        heading = ["chrom", "pos"] + stats
        outwriter.writerow(heading)
        for varfile in varfiles:
            with open(varfile, "r") as f_var:
                reader = vcf.Reader(f_var)
                for record in reader:
                    if sum([stat not in record.INFO.keys() for stat in stats]) == 0:
                        row = [record.CHROM, record.POS] + [record.INFO[stat] for stat in stats]
                        outwriter.writerow(row)

def get_specific_stats(varfile, pos, stats, chrom):
    """
    Looks through the vcf varfile until getting to the specified position
    and returns the stats for that position
    """
    import vcf
    import csv
    import sys
    with open(varfile, "r") as f:
        reader = vcf.Reader(f)
        result = None
        for record in reader:
            if int(record.POS) == pos:
                result = []
                for stat in stats:
                    if stat in list(record.INFO.keys()):
                        result.append(record.INFO[stat])
                    else:
                        result.append("NA")
                break
        if result != None:
            return( result )
        else:
            print("didn't find: " + str(chrom) + " " + str(pos))
            print("in " + varfile)
            return( ["NA"]*len(stats) )
            # sys.exit("variant not found in file")

        

def remove_alt_read_parental(infile, outfile, readcut, parent_bams, p_check = 1.0, stat_out = "trash.txt"):
    import sys
    import vcf
    import random
    """
    Remove those records where at least one of parents has greater than or
    equal to a perscribed number of ALT allele observations in the original alignment.
    Note that this refers to ALT alleles of a particular type
    """
    print("bams to be examined: " + repr(parent_bams))
    with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
        reader = vcf.Reader(f_in)
        writer = vcf.Writer(stream = f_out, template = reader)
        total_count = 0
        kept_count = 0
        for record in reader:
            if random.random() <= p_check :
                total_count += 1
                ref_base = record.REF
                pos = record.POS
                chromosome = record.CHROM
                passes_count = True
                for bam_set in parent_bams:
                    print( bam_set )
                    ## it is necessary to take pos-1 because the counting in pysam seems to start at 0?
                    record_base_counts = base_counts(bam_set, chromosome, pos-1)
                    bases = record_base_counts.keys()
                    for base in bases:
                        if str(base) != str(ref_base):
                            if record_base_counts[base] >= readcut:
                                passes_count = False
                                break
                    if not passes_count:
                        break
                if passes_count:
                    kept_count += 1
                    writer.write_record(record)
        writer.close()
    with open(stat_out, "w") as f:
        f.write(str(total_count) + "\n" + str(kept_count) + "\n")
    return (total_count, kept_count)

def base_counts(bams, chrom, position):
    """
    Count the number of times each base is found in BAMS at the specified
    POSITION on the given CHROM
    """
    import pysam
    bases = []
    if type(bams) is not list:
        bams = [bams]
    for bam in bams:
        with pysam.AlignmentFile(bam, "rb") as samfile:
            for read in samfile.fetch(chrom, start = position - 1, end = position):
                try:
                    ind = read.get_reference_positions().index(position)
                    bases.append(read.query_alignment_sequence[ind])
                except ValueError:
                    print("Warning: position not present in read, possible indel in read")
    bases_found = list(set(bases))

    result = dict()
    for base in bases_found:
        result[base] = bases.count(base)
    return result

def flanking_seq(csv_fname, outfile, ref, flank = 200, cutoff = 0.3):
    import sys
    import csv
    import pysam

    ref_fasta = pysam.FastaFile(ref)
    with open(csv_fname, "r") as f, open(outfile, "w") as fout:
        cr = csv.reader(f, delimiter = ",")
        next(cr)
        for entry in cr:
            chrom = entry[0]
            pos = int(entry[1])
            pp = float(entry[5])
            if pp > cutoff:
                fout.write('>' + chrom + '_pos_' + str(pos) + '\n')
                fout.write(ref_fasta.fetch(chrom, pos - 1 - flank, pos + flank) + '\n')


def count_gc_all(fnames, outfile, colnames, check_percs):
    import sys
    import vcf
    import random
    import gzip
    import csv

    with open(outfile, 'w') as csvfile:
        gc_writer = csv.writer(csvfile, delimiter=',')
        gc_writer.writerow(['filter', 'gc_perc'])
        for ii, fname in enumerate(fnames):
            p_check = check_percs[ii]
            reader = vcf.Reader(filename=fname, compressed=True)
            n_sites = 0.0
            n_gc = 0.0
            for record in reader:
                if random.random() <= p_check:
                    n_sites += 1
                    if record.REF == 'G' or record.REF == 'C':
                        n_gc += 1
            gc_writer.writerow([colnames[ii], n_gc/n_sites])

    

def count_dns(infile, outfile, ref, p_check = 1.0):
    import sys
    import vcf
    import random
    import gzip
    import shutil
    import csv
    import pysam

    dn_dict = dict()
    reader = vcf.Reader(filename=infile, compressed=True)
    ref_fasta = pysam.FastaFile(ref)
    for record in reader:
        if random.random() <= p_check :
            record_dn = ref_fasta.fetch(record.CHROM, record.POS - 1, record.POS + 1)
            if record_dn in dn_dict.keys():
                dn_dict[record_dn] += 1
            else:
                dn_dict[record_dn] = 1
            
    out_header = []
    for dn_type in dn_dict.keys():
        out_header.append(dn_type)
    out_counts = []
    for dn_type in out_header:
        out_counts.append(dn_dict[dn_type])
        
    with open(outfile, "w") as csvfile:
        stupid_writer = csv.writer(csvfile, delimiter = ",")
        stupid_writer.writerow(out_header)
        stupid_writer.writerow(out_counts)
        

def mask_the_vcf(infile, maskfile, outfile, stat_out = "trash.txt"):
    import sys
    import vcf
    import random
    import shutil
    import pysam

    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]

    mask_fasta = pysam.FastaFile(maskfile)
    total_count = 0
    kept_count = 0
    with open(outfile_write, 'w') as f_out:
        reader = vcf.Reader(filename=infile, compressed=True)
        writer = vcf.Writer(stream = f_out, template = reader)

        
        for record in reader:
            total_count += 1
            if (mask_fasta.fetch(str(record.CHROM), 
                                int(record.POS) - 1, 
                                int(record.POS)) == '0'):
                kept_count += 1
                writer.write_record(record)

    if os.path.splitext(outfile)[1] == '.gz':
        print("compressing output")
        with open(outfile_write, 'rb') as f_in:
            with gzip.open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(outfile_write)

    with open(stat_out, "w") as f:
        f.write(str(total_count) + "\n" + str(kept_count) + "\n")

    return(total_count, kept_count)
    

def remove_indel_parental(infile, outfile, indelcut, bams, p_check = 1.0, stat_out = "trash.txt"):
    import sys
    import vcf
    import random
    import gzip
    import shutil

    """
    Remove those variants where reads mapping to the variant sites have
    number of indel REGIONS (not sites) greater than or equal to indelcut
    """
    
    outfile_write = outfile
    if os.path.splitext(outfile)[1] == '.gz':
        outfile_write = os.path.splitext(outfile)[0]

    with open(outfile_write, 'w') as f_out:
        reader = vcf.Reader(filename=infile, compressed=True)
        writer = vcf.Writer(stream = f_out, template = reader)
        record_count = 0
        kept_count = 0
        for record in reader:
            if random.random() <= p_check :
                record_count += 1
                total_count = 0
                pos = record.POS
                chromosome = record.CHROM
                for bam_set in bams:
                    total_count += indel_count(bam_set, chromosome, pos)
                if total_count < indelcut:
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

def indel_count(bams, chrom, position):
    """
    Count the number of times an indel exists in one of the reads mapping
    to the specified site
    """
    import pysam
    if type(bams) is not list:
        bams = [bams]
    indel_count = 0
    for bam in bams:
        with pysam.AlignmentFile(bam, "rb") as samfile:
            for pileupcolumn in samfile.pileup(chrom, position-1, position):
                if pileupcolumn.pos == (position-1):
                    for pileupread in pileupcolumn.pileups:
                        cigar_stats = pileupread.alignment.get_cigar_stats()
                        ## Add up the num of regions in alignment that are insertions or deletions
                        indel_count += cigar_stats[1][1] + cigar_stats[1][2]
    return(indel_count)

def get_flank(chrom, pos, seqs, flank):
    p1 = max([pos - 1 - flank, 0])
    p2 = min([pos + flank, len(seqs[chrom])])
    return str(seqs[chrom].seq[p1:p2])

def check_repeat(seq, k):
    result = False
    sets = [[j, len(list(g))] for j, g in itertools.groupby(seq)]
    for s in sets:
        if s[1] > k and s[0] in ["A", "C", "G", "T", "a", "c", "g", "t"]:
            result = True
    return result

def count_flank_variants(chrom, pos, varfile, flank):
    positions = list(range(pos-flank, pos)) + list(range(pos+1, pos+flank+1))
    result = 0
    for position in positions:
        command = 'grep -P "' + chrom + '\t' + str(position) + '\t" ' + varfile
        args = shlex.split(command)
        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        tmp = p.communicate()[0].decode('utf-8')
        print( tmp )
        if tmp != "":
            result += 1
    return result

def count_parental_alt(record, mother, father):
    result = None
    warning = False
    m_rec = None
    for samp in record.samples:
        if samp.sample == mother:
            m_rec = samp 
    f_rec = None
    for samp in record.samples:
        if samp.sample == father:
            f_rec = samp 
    if m_rec != None or f_rec != None:
        result = m_rec.data.AD[1]
        result += f_rec.data.AD[1]
        if len(m_rec.data.AD) > 2:
            warning = True
    return result, warning

def count_alts(record, indivs):
    result = None
    warning = False
    indiv_samps = []
    for indiv in indivs:
        for samp in record.samples:
            if samp.sample == indiv:
                indiv_samps.append(samp)
    if not None in indiv_samps:
        result = 0
        for indiv_samp in indiv_samps:
            result += indiv_samp.data.AD[1]
            if len(indiv_samp.data.AD) > 2:
                warning = True
    return result, warning

def get_depths(record):
    return [samp.data.DP for samp in record.samples]

def get_parent_gt(record, mother, father):
    return [samp.data.GT for samp in record.samples if samp.sample in [mother, father]]

def get_gts(record, indivs):
    return [samp.data.GT for samp in record.samples if samp.sample in indivs]
