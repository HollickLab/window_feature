#!/usr/bin/python
'''window_features.py uses SAM file alignments to output a table of read counts in sliding windows.
Copyright (C) 2013 Joy-El R.B. Talbot
    
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
                    
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
                                    
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
	
import sys, math, array, re
import argparse		# to parse the command line arguments
import datetime     # to add dates to standard output names
import subprocess

PROGRAM = "window_feature.py"
VERSION = "0.1.1"
UPDATED = "130510 JRBT"

# RAM constants in bytes (for chunk_sizes)
ONEGB = (1024 * 1024 * 1024)
ONEMB = (1024 * 1024)
ONEKB = 1024

# Line constant = # of lines to chunk per output
# constant amounts for standard computer (16GB RAM iMac)
# update below if --highRAM option envoked
LINE = (1000 * 10)
USERAM = ONEGB

#****** Custom classes for error handling ******
class Error(Exception):
    pass

#****** Window and Step size related errors ******
class FactorError(Error):
    def __init__(self, x,y):
        self.x = x
        self.y = y
    def __str__(self):
        return "FactorError: {} is not a factor of {}.".format(self.y,self.x)

class ZeroError(Error):
    def __init__(self,name):
        self.name = name
    def __str__(self):
        return "ZeroError: A {} value of zero is not allowed.".format(self.name)

#****** End Window and Step size related errors ******
#****** End Custom classes for error handling ******

#****** Custom functions ******
def getInput(exception):
    '''Parse Error exception to return the original user input.'''

    p = str(exception).split(":")
    return p[-1].strip()


def verifyWS(w,s):
    '''Verify user inputted window (w) and step (s) values.
    
    They must be non-zero integers, and the window should be evenly 
    divisible by the step size.'''

    if w == 0 or s == 0 or w % s !=0:
        # Print a specific message to the user
        print "\nYour window size of {} and/or your step size of {} are invalid.\n".format(w,s)
        print "Both must be non-zero whole numbers and "
        print "the step size must be a factor of the window size."
        print "Please re-enter these values.\n"

        # Continue to prompt for new values until acceptable ones are given.
        while True:
            try:
                w = int(raw_input("Window Size: "))
                if w == 0: raise ZeroError("window")
                s = int(raw_input("Step Size: "))
                if s == 0: raise ZeroError("step")
                if w % s != 0: raise FactorError(w,s)
                break
            except ValueError as e:
                print "ValueError: {} is not an integer.".format(getInput(e))		
            except (FactorError, ZeroError) as e:
                print e
	
    print "The window size is: {}\nThe step size is: {}".format(w,s)
    return (w,s)


def getChromosomes(f):
    '''Confirm that SAM file is properly formatted and get the chromosome sizes.'''
    
    c = {}
    f.seek(0) # make sure you start at the top of the file
    while True:
        line = f.readline().strip()
        if len(line) == 0: break
		
        if line[0:3] == "@SQ":
            parts = line.split("\t")
            c[parts[1].strip()[3:]] = int(parts[2].strip()[3:])
		
        if line[0:1] != "@":
            # finished with header, so do some checks
            parts = line.split("\t")
            # have @SQ header lines?
            if len(c.keys()) == 0:
                sys.exit("No chromosome sizes provided in @SQ header lines.")
            # have all required SAM columns?
            if len(parts) < 11:
                sys.exit("Incomplete SAM file, missing required columns.")
            # if using -a or -m options do you have the appropriate option columns?
            tags = {}
            for p in range(11,len(parts)):
                tags[parts[p][0:2]] = 0
            if args.normalize:
                if "NH" not in tags:
                    sys.exit("Missing NH:i:00 tags required for the --multireads option.")
            if args.abundance:
                if "NA" not in tags:
                    sys.exit("Missing NA:i:00 tags required for the --abundance option.")
            break
    f.seek(0) # restart file at first line
    return c


def read_chunk(file_obj,chunk_size):
    '''Read in file by chunk_size chunks returning one line at a time.'''
    chunk_num = 1
    sys.stdout.write(".")
    sys.stdout.flush()
    # get first chunk
    chunk = file_obj.read(chunk_size)
    # continue looping until a chunk is just EOF (empty line)
    while chunk:
        chunk_num += 1
        chunk_list = chunk.split("\n")
        # yield all but last, potentially incomplete line
        for c in chunk_list[:-1]:
            yield c
        sys.stdout.write(".")
        sys.stdout.flush()
        if (chunk_num % 100) == 0:
            sys.stdout.write(" {}\n".format(chunk_num))
            sys.stdout.flush()
        # add incomplete line to beginning of next chunk read
        chunk = chunk_list[-1] + file_obj.read(chunk_size)
    sys.stdout.write("\n")
    sys.stdout.flush()

def getLength(cigar):
    '''Parse out the length of the read from the CIGAR notation of its alignment.
    
    Number values from 00M, 00X, 00I, 00S, 00= should all be added together to get the total length.'''
    
    number = ""
    count = 0
    for j in range(len(cigar)):
        try:
            
            number += str(int(cigar[j]))
        except ValueError:
            if cigar[j] == 'M' or cigar[j] == 'X' or cigar[j] == 'I' or cigar[j] == 'S' or cigar[j] == '=':
                count += int(number)
            number = 0
    return count

def get_arguments(program, version, update):
    date = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    
    parser = argparse.ArgumentParser(description='''Count reads per window and output results to a tab-delimited file.''')

    parser.add_argument("-v","--version", 
	    				action  = 'version', 
		    			version = "{} {}\tUpdated {}".format(program, version, update))
    parser.add_argument("--sam",
                        help    = "input alignment file in SAM format",
                        metavar = 'FILE')
    parser.add_argument("--bam",
                        help    = "input alignment file in BAM format",
                        metavar = 'FILE')
    parser.add_argument("-o", "--output", 
	    				help    = "output file name",
		    			default = date + "_windows.csv")
    parser.add_argument("-w","--window", 
	    				help    = "define window-size in bases (default = 1,000,000)", 
		    			default = 1000000, 
			    		metavar = 'N',
				    	type    = int)
    parser.add_argument("-s","--step", 
	    				help    = "define step-size for sliding window in bases (default = 10,000)", 
		    			default = 10000, 
			    		metavar = 'N',
				    	type    = int)
    parser.add_argument("-n","--normalize", 
	    				help    = "normalize multi-reads as fractional hits per alignment", 
		    			action  = 'store_true')
    parser.add_argument("-a","--abundance", 
	    				help    = "calculate read abundance from SAM options", 
		    			action  = 'store_true')
    parser.add_argument("--sizemin", 
                        help    = "Minimum size read to tally, default = 20",
                        default = 20,
                        type    = int,
                        metavar = 'N')
    parser.add_argument("--sizemax",
                        help    = "Maximum size read to tally, default = 25",
                        default = 25,
                        type    = int,
                        metavar = 'N')
    parser.add_argument("--poolstrands",
                        help    = "Pool tallies on both strands, default = false",
                        action  = 'store_true')
    
    return parser.parse_args()

def header(filename):
    '''Read in header with samtools view -H and return chromosome info

    Later use to create chromosome objects, for now just return a dictionary of chromosome names + lengths.'''

    if filename[-3:] == "sam":
        headerlines = subprocess.check_output(["samtools","view","-H","-S",filename], stderr=subprocess.PIPE).split("\n") # stderr arg hides stderr from user
    elif filename[-3:] == "bam":
        headerlines = subprocess.check_output(["samtools","view","-H",filename], stderr=subprocess.PIPE).split("\n")
    else:
        sys.exit("{} is not a SAM or BAM file".format(filename))

    # Later implement function checking using the stderr PIPE output
    
    chrs = {}

    for l in headerlines:
        if l[0:3] == "@SQ":
            try:
                name = re.search('SN:(\S+)',l).group(1)
                length = int(re.search('LN:(\d+)',l).group(1))
                yield (name,length)
            except:
                raise Error(l, "Could not extract name and length data from this SAM header:\t{}".format(l))


def bam_chunk(filename, chunksize, chromosomes):
    '''Parse through a bam file in chunks and return line-by-line. 
    
    Importantly, chunks are measured in nucleotides!'''
    
    checkIndexBam(filename)
    
    for c in sorted(chromosomes.keys()):
        start = 1 # because sam file is in 1-based

        while start < int(chromosomes[c]):
            # loop through entire chromosome
            # Status tracker 
            sys.stdout.write(".")
            if ((start + chunksize) / chunksize) % 100 == 0:
                sys.stdout.write("{} chunks".format(start / chunksize))
            sys.stdout.flush()
            #print "{}:{}-{} group:".format(c, start, start + chunksize - 1)
            lines = subprocess.check_output(["samtools", "view", filename, "{}:{}-{}".format(c, start, start + chunksize - 1)]).split("\n")
            for l in lines:
                if len(l) > 0:
                    sam = l.split("\t")
                    #sys.stdout.write(sam[0])
                    # only return lines that START in the current chunk
                    if int(sam[3]) >= start:
                        #sys.stdout.write(": PASSES - {}".format(sam[3]))
                        yield sam
                    #else: sys.stdout.write(": FAILED - {}\n".format(sam[3]))
                    #sys.stdout.flush()
            # loop start
            start += chunksize
        sys.stdout.write("\n")
        sys.stdout.flush()

def checkIndexBam(filename):
    '''Confirm that file is in the bam format and has an index file.'''
    
    try:
        subprocess.check_output(["samtools", "view", "-H", filename], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        sys.exit("{} is not in a recognizable bam format.".format(filename))

    try:
        subprocess.check_output(["samtools", "view", filename, "Chr1:100-200"], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        sys.exit("Could not find an index (.bai) file for {}.".format(filename))

    return True



#****** Validate user input ******
starttime = datetime.datetime.now()
print "Start program: {}".format(starttime.strftime('%H:%M:%S'))

args = get_arguments(PROGRAM, VERSION, UPDATED)

# Make sure the input file is a valid file
if args.sam:
    infile = args.sam
elif args.bam:
    infile = args.bam
else:
    sys.exit("You must include an input file (SAM or BAM).")

try:
    inf = open(infile)
except IOError as e:
    sys.exit("{} is not an existing file or directory.".format(getInput(e)))

# Make sure that output file doesn't exist
try:
    outf = open(args.output)
except IOError as e:
    if e.strerror == "No such file or directory":
        # only create the file if a "No such file or directory" IOError occurred
        outf = open(args.output, 'w')
    else:
        raise IOError(e)
else:
    sys.exit("{} output file already exists.".format(args.output))

# define dictionary of chromosome keys and length values from SAM file
print "Load chromosome info from {}: {}".format(infile, datetime.datetime.now().strftime('%H:%M:%S'))
chromosomes = {}
for c in header(infile):  # works for both SAM and BAM files
    chromosomes[c[0]] = c[1]

# Make sure step is a factor of window and both are non-zero integers
(window, step) = verifyWS(args.window, args.step)

# Now to start doing the work
# array of nucleotide sizes to keep
read_size = {}
for r in range(args.sizemin,args.sizemax + 1):
    read_size[r] = 1

# initialize steps
print "Initializing step bins: {}".format(datetime.datetime.now().strftime('%H:%M:%S'))
steps = {}

for c in chromosomes.keys():
    if args.poolstrands: strand = 1
    else: strand = 2
    
    steps[c] = [[[0 for s in range(strand)] for l in range(len(read_size))] for b in range(int(math.ceil(chromosomes[c] / float(step))))]
        
# start looping through the file
print "Starting to tally alignments into bins: {}".format(datetime.datetime.now().strftime('%H:%M:%S'))
print "Each . represents {} lines of the samfile being read and processed.".format(ONEMB)

if args.bam:
    for parts in bam_chunk(infile, ONEMB, chromosomes):
        # add sam line to window counts
        r_length = getLength(parts[5])
        if r_length in read_size:
            
            r_chromosome = parts[2]
            # determine strandedness
            strand = 0
            if not args.poolstrands and int(parts[1]) == 16:
                strand = 1 
        
            r_start = int(parts[3])
            r_end = r_start + r_length - 1
            
            # deal with abundance and normalization
            if args.abundance:
                for i in range(len(parts)-1,-1,-1): # loop backwards
                    if parts[i][0:2] == "NA":
                        na_parts = parts[i].split(":")
                        na = int(na_parts[-1])
                        break
            else:
                na = 1
            
            if args.normalize:
                for i in range(len(parts)-1,-1,-1):
                    if parts[i][0:2] == "NH":
                        nh_parts = parts[i].split(":")
                        nh = float(nh_parts[-1])
                        break
            else: 
                nh = 1
             
            # must make step a float in order for the math to work
            bin_s = int(math.ceil(r_start / float(step))) - 1
            bin_e = int(math.ceil(r_end / float(step))) - 1
        
            if bin_s == bin_e:
                steps[r_chromosome][bin_s][(r_length - args.sizemin)][strand] += ((r_length * na) / nh)
            else:
                # add length in bin_s
                steps[r_chromosome][bin_s][(r_length - args.sizemin)][strand] += (((int((bin_s + 1) * step) - r_start + 1) * na) / nh)
                steps[r_chromosome][bin_e][(r_length - args.sizemin)][strand] += (((r_end - int(bin_e * step)) * na ) / nh)
            
                # if more than two bins are spanned
                full_bins = bin_e - bin_s - 1

                for b in range(1,full_bins + 1):
                    steps[r_chromosome][bin_s + b][(r_length - args.sizemin)][strand] += ((step * na) / nh) 
      
elif args.sam:
    for line in read_chunk(inf, ONEMB):
        line = line.strip()
                
        parts = line.split("\t")
        
        # only process non-comment, mapping alignments
        if len(line) != 0 and line[0] != "@" and int(parts[1]) != 4:
            # add sam line to window counts
            r_length = getLength(parts[5])
            if r_length in read_size:
            
                r_chromosome = parts[2]
                # determine strandedness
                strand = 0
                if not args.poolstrands and int(parts[1]) == 16:
                    strand = 1 
        
                r_start = int(parts[3])
                r_end = r_start + r_length - 1
            
                # deal with abundance and normalization
                if args.abundance:
                    for i in range(len(parts)-1,-1,-1): # loop backwards
                        if parts[i][0:2] == "NA":
                            na_parts = parts[i].split(":")
                            na = int(na_parts[-1])
                            break
                else:
                    na = 1
            
                if args.normalize:
                    for i in range(len(parts)-1,-1,-1):
                        if parts[i][0:2] == "NH":
                            nh_parts = parts[i].split(":")
                            nh = float(nh_parts[-1])
                            break
                else: 
                    nh = 1
             
                # must make step a float in order for the math to work
                bin_s = int(math.ceil(r_start / float(step))) - 1
                bin_e = int(math.ceil(r_end / float(step))) - 1
        
                if bin_s == bin_e:
                    steps[r_chromosome][bin_s][(r_length - args.sizemin)][strand] += ((r_length * na) / nh)
                else:
                    # add length in bin_s
                    steps[r_chromosome][bin_s][(r_length - args.sizemin)][strand] += (((int((bin_s + 1) * step) - r_start + 1) * na) / nh)
                    steps[r_chromosome][bin_e][(r_length - args.sizemin)][strand] += (((r_end - int(bin_e * step)) * na ) / nh)
            
                    # if more than two bins are spanned
                    full_bins = bin_e - bin_s - 1

                    for b in range(1,full_bins + 1):
                        steps[r_chromosome][bin_s + b][(r_length - args.sizemin)][strand] += ((step * na) / nh) 

inf.close()

# define windows and group bins into windows
print "Now combining bins into Windows: {}".format(datetime.datetime.now().strftime('%H:%M:%S'))
header = "Chromosome\tWindow"
if args.poolstrands:
    for d in range(args.sizemin,args.sizemax +1):
        header += "\t" + str(d)
else:
    for d in range(args.sizemin, args.sizemax +1):
        header += "\t" + str(d) + "W" + "\t" + str(d) + "C"
outf.write(header + "\n")

windows = {}
window_count = int(math.ceil((chromosomes[c] - window) / float(step))) + 1
for c in chromosomes.keys():
    if args.poolstrands: strand = 1
    else: strand = 2
    
    windows[c] = [[[0 for s in range(strand)] for l in range(len(read_size))] for b in range(window_count)]
    
    for w in range(window_count):
        for s in range(int(window / step)):
            for l in range(len(read_size)):
                for x in range(strand):
                    # w+s gives bin to use, w = # bins before we start counting, s = bin in current window 
                    windows[c][w][l][x] += steps[c][w+s][l][x]

    for w in range(window_count):
        output = str(c) + "\t" + str(w)
        for l in range(len(read_size)):
            for x in range(strand):
                output += "\t" + str(windows[c][w][l][x])
        outf.write(output + "\n")

outf.close()
inf.close()
endtime = datetime.datetime.now()
elapsed = endtime - starttime
e_seconds = elapsed.total_seconds()
e_hours = int(e_seconds / 3600)
e_minutes = int( (e_seconds - (e_hours * 3600)) / 60)
e_seconds = (e_seconds - (e_hours * 3600) - (e_minutes * 60))
print "{} run finished: {}".format(PROGRAM, endtime.strftime('%H:%M:%S'))
print "Total runtime: {} hours {} minutes {} seconds".format(e_hours, e_minutes, e_seconds)
