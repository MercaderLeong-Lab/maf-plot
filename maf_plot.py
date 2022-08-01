import argparse
import pandas as pd
import seaborn as sns
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt

#The default expected header of a snp-stats input file
HEADER = "alternate_ids rsid chromosome position alleleA alleleB comment HW_exact_p_value HW_lrt_p_value alleleA_count alleleB_count alleleA_frequency alleleB_frequency minor_allele_frequency minor_allele major_allele info impute_info missing_proportion A B AA AB BB NULL total"

#A alternative header for input files without Hardy-Weinberg values
ALT_HEADER="alternate_ids rsid chromosome position alleleA alleleB comment alleleA_count alleleB_count alleleA_frequency alleleB_frequency minor_allele_frequency minor_allele major_allele info impute_info missing_proportion A B AA AB BB NULL total"

def use_alt_header():
    '''Switch to using the alternate header'''
    global COLUMN, HEADER
    HEADER = ALT_HEADER
    COLUMN = {value: location for location, value in enumerate(ALT_HEADER.split())}
    print("Using Alternate Header")

#A map from SNP value to the column it is located in
COLUMN = {value: location for location, value in enumerate(HEADER.split())}

class MAF_Reader:
    '''A wrapper class for a file that contains MAF data.
       The file should not have a header and the data should be organized as follows: chromosome, position, maf, major allele, minor allele.
       The file should be sorted first by numericla Chromosome number then by position within chromosome.
    '''

    def __init__(self, source):
        self.source = source
        self.chr = None 
        self.ended = False

    def advance(self):
        '''Advances the reader a single line and parses the new line. The values are saved as instance variables'''
        line = self.source.readline().split()
        if line == "":
            self.ended = True
            return
        self.chr = line[0]
        self.pos = int(line[1])
        self.maf = float(line[2])
        self.maj = line[3]
        self.min = line[4]

    def advance_to(self, chr):
        '''Advances the reader to the given chromosome. No guarantees about what happens if the chromosome is before the current one or if the input is an invalid chromosome name.'''
        while self.chr != chr and not self.ended:
            self.advance()
        
class Entry:
    '''A class that represents a single line entry of a snp stats file'''

    def __init__(self, data):
        '''Parses a line of text into a snp stats line. Optionally, a MAF_Source object may be provided.'''
        self.data = data

    def to_string(self):
        return " ".join(self.data) + "\n"

    def info(self):
        '''Returns the info value of the variant or -1 if it cannot be parsed.'''
        try:
            return float(self.data[COLUMN["info"]])
        except ValueError:
            return -1

    def chr(self):
        '''Returns the chromosome that the variants is located on'''
        return self.data[COLUMN["chromosome"]]

    def pos(self):
        '''Returns the positions within the chromosome of the variant'''
        return int(self.data[COLUMN['position']])

    def maf(self):
        '''Returns the MAF of the variants, either from the snp stats line or a separate file if provided.'''
        try:
            return float(self.data[COLUMN["minor_allele_frequency"]])
        except ValueError:
            return -1

    def set_maf(self, maf):
        self.data[COLUMN["minor_allele_frequency"]] = str(maf)

    def missing_proportion(self):
        return float(self.data[COLUMN["missing_proportion"]])

    def structural_variant(self):
        '''Returns true if the variant is a structural variant. A structural variants is defined as any variant where either allele is not a single nucleotide'''
        return len(self.data[COLUMN["alleleA"]]) != 1 or len(self.data[COLUMN["alleleB"]]) != 1

    @staticmethod
    def for_entry(source, f):
        '''Iterates over an entire snp stats file and runs f on each data line of the file'''
        entries = 0
        for line in source:
            if line[0] == "#": continue
            if HEADER in line: continue
            if ALT_HEADER in line: 
                use_alt_header()
                continue
            f(Entry(line.split()))

    @staticmethod
    def for_folder(folder, f, g = lambda x: None):
        '''Iterates over an entire folder an runs f on each line in a snp stats file in the folder. The folder must have files for all 22 autosomes labeled chrN-snp-stats.txt. The function f is run for each entry in each file. The function g is run once before each chromosome file is read.'''
        for i in range(1, 23):
            g(i)
            Entry.for_entry(open(folder + '/' + file_name(i), 'r'), f)

class Bin:
    '''A utility class that represents a 'bin' of items within a range where the average and count are important but each individual value is not essential.'''
    def __init__(self, minimum, maximum):
        self.minimum = minimum
        self.maximum = maximum
        self.n = 0
        self.average = 0

    def __str__(self):
        return str(self.minimum) + "-" + str(self.maximum)

    def add(self, value):
        '''Adds the given value to the bin, updating the count and average.'''
        self.n += 1
        self.average = self.average + (value - self.average)/self.n

    def in_range(self, item):
        '''Returns whether the given value is in the range of the bin'''
        return (self.minimum < item) and (item < self.maximum)

    def add_if_in_range(self, x, y):
        '''Adds y is x is in the range of the bin and returns whether the values was added.'''
        if self.in_range(x):
            self.add(y)
            return True
        else:
            return False

    def average(self):
        '''Returns the current average of the bin'''
        return self.average

    def count(self):
        '''Returns the number of items in the bin'''
        return self.n

    @staticmethod
    def default_bins():
        '''Returns a list of default bins that cover the range [0, 0.5]'''
        return [
            Bin(0, 0.0005),
            Bin(0.0005, 0.001),
            Bin(0.001, 0.002),
            Bin(0.002, 0.005),
            Bin(0.005, 0.01),
            Bin(0.01, 0.015),
            Bin(0.015, 0.02),
            Bin(0.02, 0.035),
            Bin(0.035, 0.05),
            Bin(0.05, 0.1),
            Bin(0.1, 0.2),
            Bin(0.2, 0.3),
            Bin(0.3, 0.4),
            Bin(0.4, 0.5)
        ]

    @staticmethod
    def find_bin(bins, item):
        '''Returns the first bin that the items belongs in.'''
        for bin in bins:
            if bin.in_range(item): return bin

    @staticmethod
    def sort(bins, key, item):
        '''Inserts item into the first bin that key is in the range of. Returns whether the item was added to a bin.'''
        for bin in bins:
            if bin.add_if_in_range(key, item): return True
        return False

def write_header(f):
    f.write(HEADER + "\n")

def file_name(chrom, prefix="", suffix="-snp-stats.txt"):
    '''Returns the expected snp-stats file name for a file'''
    return "chr" + str(chrom) + "-snp-stats.txt"

hits = 0
misses = 0
dest_file = None

def remaf(source, maffile, dest):
    '''Copies files from the source folder to the destination folder, replace variant MAF with the MAF in the given file or -1 if not found. The format of the MAF file must be CHR POS MAF MAJ MIN'''
    global hits
    global misses
    global dest_file
    print("Running Remaf")

    print("Loading MAF Reader")
    maf_reader = MAF_Reader(maffile)
    print("MAF Reader Loaded")

    def process_entry(entry):
        global hits
        global misses
        global dest_file
        pos = entry.pos()
        chr = entry.chr()

        while not maf_reader.ended and maf_reader.chr == chr and maf_reader.pos < pos:
            maf_reader.advance()

        if not maf_reader.ended and maf_reader.chr == chr and maf_reader.pos == pos:
            hits += 1
            entry.set_maf(maf_reader.maf)
        else:
            misses += 1
            entry.set_maf(-1)

        dest_file.write(entry.to_string())

    for i in range(1,23):
        print(f"Loading chr{i}")
        maf_reader.advance_to("chr" + str(i))
        src_file = open(source + "/" + file_name(i), 'r')
        dest_file = open(dest + "/" + file_name(i), 'w')

        write_header(dest_file)

        hits = 0
        misses = 0

        Entry.for_entry(src_file, process_entry)

        print(f"Loaded with {hits} hits and {misses} misses")

        src_file.close()
        dest_file.close()

def lineplot(source, output):
    '''Generates a lineplot with a line for each source folder comparing MAF and info score'''
    print("Generating Lineplot")

    bins = Bin.default_bins()

    stats = {'info': [], 'bin': [], 'panel': []}

    for folder in args.source:

        def process_entry(entry):
            nonlocal stats, bins
            info = entry.info()
            maf = entry.maf()
            if info == -1: return
            if maf == -1: return
            stats['info'].append(info)
            stats['bin'].append(str(Bin.find_bin(bins, maf)))
            stats['panel'].append(folder)
        
        print(f'Loading {folder}')
        Entry.for_folder(folder, process_entry, lambda chr: print(f"Reading chr{chr}"))

    stats = pd.DataFrame(stats)

    order = CategoricalDtype(
        [str(b) for b in bins],
        ordered=True
    )

    stats['bin'] = stats['bin'].astype(order)

    print("Generating Plot")

    sns.set_style("whitegrid")

    plot = sns.lineplot(data = stats.sort_values('bin'), x='bin', y='info', hue='panel', sort=False)

    plot.set(xlabel="MAF", ylabel="Info")
    plot.set(title="Average Info by MAF")
    plt.xticks(rotation=75)
    plt.grid(b=True)
    plt.tight_layout()

    plt.savefig(args.output, format="svg")

    print("Done")

def barplot(source, output, thresholds, y_range=None):
    print("Generating Barplot")
        
    n = len(args.thresholds)
    ranges = [None] * n
    ranges[0] = vc.Bin(0, args.thresholds[0])
    
    for i in range(len(args.thresholds) - 1): ranges[i + 1] = vc.Bin(args.thresholds[i], args.thresholds[i + 1])

    ranges.append(vc.Bin(args.thresholds[-1], 1))
    
    bins = Bin.default_bins()
    categories = [str(b): vc.Bin.default_bins() for b in ranges } 

    def update_count(entry):
        nonlocal categories
        info = entry.info()
        maf = entry.maf()
        if info == -1: return
        if maf == -1: return
        for band in ranges:
            if band.in_range(enty.info()):
                Bin.sort(categories[str(band)], maf, info)

    Entry.for_folder(source, update_count, lambda x: print(f"Reading chr{x}"))

    print("Generating Plot")

    bars = pd.DataFrame({band: [b.count() for b in bins] for (band, bins) in categories.items{}}, index = [str(b) for b in bins])

parser = argparse.ArgumentParser()
parser.add_argument('-p', nargs='?', default="", type=str, help="The file prefix before the chromosome number of each snp-stats file")
parser.add_argument('-P', nargs='?', default="-snp-stats.txt", type=str, help="The file suffix bfore the chromosome number of each snp-stats file")
parser.add_argument('-r', nargs=2, default=[1,22], type=int)

subparsers = parser.add_subparsers(dest="command", required=True)

remaf_parser = subparsers.add_parser("remaf")
lineplot_parser = subparsers.add_parser("lineplot")
barplot_parser = subparsers.add_parser("barplot")
test_parser = subparsers.add_parser("test")

remaf_parser.add_argument("source", type=str)
remaf_parser.add_argument("maffile", type=argparse.FileType('r'))
remaf_parser.add_argument("dest", type=str)

lineplot_parser.add_argument("source", type=str, nargs="+", help="The folders containing the snp-stats data")
lineplot_parser.add_argument("output", type=argparse.FileType('w'), help="The file to write the graph to")

barplot_parser.add_argument("source", type=str, help="The folder containing snp-stats data")
barplot_parser.add_argument("output", type=argparser.FileType('w'), help="The image file to save the graph in")
barplot_parser.add_argument("thresholds", type=float, nargs="+", help="The info thresholds to divide the variants by")


if __name__ == "__main__":
    args = parser.parse_args()
    chromosomes = range(args.r[0], args.r[1] + 1)

    if args.command == "remaf":
        remaf(args.source, args.maffile, args.dest)

    elif args.command == "lineplot":
        lineplot(args.source, args.output)

    elif args.command == "barplot":
        barplot(args.source, args.output, args.thresholds, y_range= args.m if args.m else None)

    elif args.command == "test":
        print(args)
        print("Testing done")
