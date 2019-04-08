import csv as cessvee
import os
import sys
import time
import subprocess
import shlex
import signal
from platform import python_version

def make_bam_plots(bams, pos, chrom, outdir, genome,
                   region=200, collapse = True,
                   height=3000, width=1600, max_panel_height=1000):

    print("Working in : " + outdir + "\n")
    sys.stderr.write("Working in : " + outdir)
    command_base = "igv_plotter"
    command_base += " --height " + str(height) + " --width " + str(width)
    command_base += " --max-panel-height " + str(max_panel_height)
    command_base += " -o " + outdir
    command_base += " -v"
    command_base += " -m 2G"
    if collapse:
        command_base += " --collapse "
    command_base += " -g " + genome
    for bam in bams:
        command_base += " " + bam
    print( "Base Plotting command to be run:" )
    print( command_base )
    complete = False
    counter = 1
    small_fname = os.path.join(outdir, "s1__" + chrom + "_" + str(pos) + ".png")
    while not complete:
        small_proc = subprocess.Popen(shlex.split(command_base + " " + chrom + ":" + str(pos)),
                                      preexec_fn=os.setsid)
        print("Trying to make small plot for time: " + str(counter))
        try:
            small_proc.wait(timeout=300)
            print("Checking that " + small_fname + " was created...")
            complete = os.path.isfile(small_fname)
            if not complete:
                try:
                    print("Attempting to kill process")
                    os.killpg(os.getpgid(small_proc.pid), signal.SIGTERM)
                except ProcessLookupError:
                    print("Process already dead")
        except subprocess.TimeoutExpired:
            print("Took too long to make figure -- killing proces...")
            try:
                print("Attempting to kill process")
                os.killpg(os.getpgid(small_proc.pid), signal.SIGTERM)
            except ProcessLookupError:
                print("Process already dead")
            counter += 1
    complete = False
    counter = 1
    big_fname = os.path.join(outdir, "s1__" + chrom + "_" + str(pos - region) + "_" + str(pos + region) + ".png")
    while not complete:
        big_proc = subprocess.Popen(shlex.split(command_base + " " + chrom + ":" +
                                                  str(pos - region) + "-" + str(pos + region)),
                                    preexec_fn=os.setsid)
        print("Trying to make big plot for time: " + str(counter))
        try:
            big_proc.wait(timeout=300)
            print("Checking that " + big_fname + " was created...")
            complete = os.path.isfile(big_fname)
            if not complete:
                try:
                    print("Process already dead")
                    os.killpg(os.getpgid(big_proc.pid), signal.SIGTERM)
                except ProcessLookupError:
                    print("Process already dead")
        except subprocess.TimeoutExpired:
            print("Took too long to make figure -- killing proces...")
            try:
                print("Attempting to kill process")
                os.killpg(os.getpgid(big_proc.pid), signal.SIGTERM)
            except ProcessLookupError:
                print("Process already dead")
            counter += 1
    print("...finished attempting to make plots\n")

def monitor_process(command, badstrings=[None]):
    print( command )
    process = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE)
    while True:
        output = process.stderr.readline().strip().decode("ascii")
        print( output )
        if True in [badstring in output for badstring in badstrings]:
            process.kill()
            sys.stderr.write( "command: " + command  + "\n")
            sys.stderr.write( "offending line: " + output + "\n")
            sys.exit("encountered a strange error!")
        if process.poll() is not None:
            break

def get_bam_region(fname, outname, chrom, pos, width=100):
    command = "samtools view -b " + fname + " " + chrom + ":" + repr(int(pos-width)) + "-" + repr(int(pos+width))
    command += " > " + outname
    os.system(command)
    os.system("samtools index " + outname)
    return command

def merge_bams(outname, bams):
    if len(bams)==1 :
        os.system("mv " + bams[0] + " " + outname)
        os.system("samtools index " + outname)
        command = "Message: no need to merge, single bam"
    else:
        command = "samtools merge -f " + outname 
        for bam in bams:
            command += " " + bam 
        os.system(command)
        for bam in bams:
            os.system("rm " + bam)
        os.system("samtools index " + outname)
    return command

def get_identifier(output_name, cname):
    identifier = os.path.basename(output_name).split("trio")[0]
    if identifier == "":
        identifier = "none"
    elif identifier[len(identifier) - 1] == "_" :
        identifier = identifier[0:(len(identifier) - 1)]
    identifier += "_" + cname
    return identifier

def get_bams_and_plot(output, csvs, bams_child, bams_mother, bams_father,
                      cname, mname, fname,
                      cutoff, outdir, identifier, region_size, genome):
    outdirs = []
    cwd = None

    for csv in csvs:
        print( csv )
        with open(csv, "r") as f:
            cr = cessvee.reader(f, delimiter=",")
            next(cr)  # skip the header
            for entry in cr:
                chrom = entry[0]
                pos = int(entry[1])
                pp = float(entry[5])
                if pp > cutoff :
                    print( entry )
                    current_outnames = []
                    ## Get child subsets ##
                    outname = (outdir + "trio_" + cname + "_" + cname + "_" + chrom +
                               "_" + repr(pos) + "_" + identifier +".bam")
                    tmpnames = []
                    counter = 1
                    # get the region in each child bam then merge them
                    for bam in bams_child:
                        tmpname = cname + "_tmp" + str(counter) + "_" + identifier + ".bam"
                        foo = get_bam_region(fname=bam, outname=tmpname, chrom=chrom, pos=pos,
                                             width=int(region_size/2))
                        print( foo )
                        tmpnames.append(tmpname)
                        counter += 1
                    foo = merge_bams(outname, tmpnames)
                    print( foo )
                    current_outnames.append(outname)
                    ## Get mother subsets ##
                    outname = (outdir + "trio_" + cname + "_" + mname + "_" + chrom +
                               "_" + repr(pos)+ "_" + identifier + ".bam")
                    tmpnames = []
                    counter = 1
                    # get the region in each child bam then merge them
                    for bam in bams_mother:
                        tmpname = cname + "_tmp" + str(counter) + "_" + identifier+ ".bam"
                        foo = get_bam_region(fname=bam, outname=tmpname, chrom=chrom, pos=pos,
                                             width=int(region_size/2))
                        print( foo )
                        tmpnames.append(tmpname)
                        counter += 1
                    foo = merge_bams(outname, tmpnames)
                    print( foo )
                    current_outnames.append(outname)
                    ## Get father subsets ##
                    outname = (outdir + "trio_" + cname + "_" + fname + "_" + chrom +
                               "_" + repr(pos)+ "_" + identifier + ".bam")
                    tmpnames = []
                    counter = 1
                    # get the region in each child bam then merge them
                    for bam in bams_father:
                        tmpname = cname + "_tmp" + str(counter)+ "_" + identifier + ".bam"
                        foo = get_bam_region(fname=bam, outname=tmpname, chrom=chrom, pos=pos,
                                             width=int(region_size/2))
                        print( foo )
                        tmpnames.append(tmpname)
                        counter += 1
                    foo = merge_bams(outname, tmpnames)
                    print( foo )
                    current_outnames.append(outname)
                    dirname = chrom + "_" + str(pos) + "_" + cname + "_" + identifier + "/"
                    outdirs.append(dirname)
                    os.chdir(outdir)
                    # test that changing directories worked
                    cwd = os.getcwd()
                    print( "cwd: " + cwd )
                    os.system("rm -rf " + dirname + " ; mkdir " + dirname)
                    command = "mv "
                    for name in current_outnames:
                        command += name + " "
                        command += name + ".bai "
                    command += dirname
                    os.system(command)
                    print( command )
                    ## MAKE THE BAM PLOTS, BECAUSE THESE ARE IN THE DIRECTORY 'DIRNAME'
                    ## THEY WILL BE INCLUDED IN TAR
                    bam_names = [dirname + os.path.basename(name) for name in current_outnames]
                    make_bam_plots(bams = bam_names,
                                   pos = pos,
                                   chrom = chrom,
                                   outdir = dirname,
                                   region = 200,
                                   genome = genome)
                    ## Write the pp of the variant to a file in this directory
                    with open(dirname + "pp_val.txt", "w") as pp_file:
                        pp_file.write(str(pp) + "\n")

    if cwd != None :
        print( "cwd: " + cwd )
        command = "tar -czf " + output + " -C " + outdir
        for d in outdirs:
            command += " " + d
            os.system(command)

        # delete the files now that the archive has been created
        for d in outdirs:
            os.system("rm -rf " + d)
    else:
        os.system("touch " + output)
    return outdirs

def get_bams_and_plot_alt(output, sitefile_name, cname, mname, fname, cutoff, bamdir, genome,
                          identifier="tmp", region_size=200):
    '''
    The sitefile must have the format:
    chromosome,position,pp
    chr01,2000,0.1
    '''
    sitefile = open(sitefile_name, "r")
    sitereader = cessvee.reader(sitefile, delimiter=",")

    outdir = bamdir

    outdirs = []
    cwd = None

    for site in sitereader:
        print( site )
        chrom = site[0]
        pos = int(site[1])
        pp = float(site[2])
        if(pp >= cutoff):
            big_bam_child = os.path.join(bamdir, chrom + "_" + cname + ".bam")
            big_bam_mother = os.path.join(bamdir, chrom + "_" + mname + ".bam")
            big_bam_father = os.path.join(bamdir, chrom + "_" + fname + ".bam")

            c_outname = (outdir + "trio_" + cname + "_" + cname + "_" + chrom +
                         "_" + repr(pos) + "_" + identifier +".bam")
            small_bam_child = get_bam_region(fname=big_bam_child,
                                             outname=c_outname, chrom=chrom, pos=pos,
                                             width=int(region_size/2))
            m_outname = (outdir + "trio_" + mname + "_" + mname + "_" + chrom +
                         "_" + repr(pos) + "_" + identifier +".bam")
            small_bam_child = get_bam_region(fname=big_bam_mother,
                                             outname=m_outname, chrom=chrom, pos=pos,
                                             width=int(region_size/2))
            f_outname = (outdir + "trio_" + fname + "_" + fname + "_" + chrom +
                         "_" + repr(pos) + "_" + identifier +".bam")
            small_bam_father = get_bam_region(fname=big_bam_father,
                                              outname=f_outname, chrom=chrom, pos=pos,
                                              width=int(region_size/2))

            ## NOW MOVE ALL THE SMALL BAM FILES TO A NEW DIRECTORY
            dirname = chrom + "_" + str(pos) + "_" + cname + "_" + identifier + "/"
            outdirs.append(dirname)
            os.chdir(outdir)
            # test that changing directories worked
            cwd = os.getcwd()
            print( "cwd: " + cwd )
            os.system("rm -rf " + dirname + " ; mkdir " + dirname)
            command = "mv " 
            for name in [c_outname, m_outname, f_outname]:
                command += name + " "
                command += name + ".bai "
                command += dirname
            os.system(command)
            print( command )

            ## MAKE THE BAM PLOTS, BECAUSE THESE ARE IN THE DIRECTORY 'DIRNAME'
            ## THEY WILL BE INCLUDED IN TAR
            bam_names = [dirname + os.path.basename(name) for name in [c_outname, m_outname, f_outname]]
            make_bam_plots(bams = bam_names,
                           pos = pos,
                           chrom = chrom,
                           outdir = dirname,
                           region = 200,
                           genome = genome)
            ## Write the pp of the variant to a file in this directory
            with open(dirname + "pp_val.txt", "w") as pp_file:
                pp_file.write(str(pp) + "\n")

            command = "tar -czf " + output + " -C " + outdir
            for d in outdirs:
                command += " " + d
            os.system(command)

            # delete the files now that the archive has been created
            for d in outdirs:
                os.system("rm -rf " + d)

        return outdirs
