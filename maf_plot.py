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

class MAF_Source:
    '''This is a class that provdies the ability to load a MAF from an external file. In order for it to work, the file must be sorted first by chromosome, then by position within the chromosome. The first three columns of the file must be chromosome number, position, and minor allele frequency.'''
    
    header = [
        'chr',
        'position',
        'maf',
    ]

    def __init__(self, source):
        '''Initializes a source file for later use. This can take some time as it caches the start position of each line in the file for quick later access.'''
        self.source = source
        self.line_offsets = [] 
        self.chr_lines = {}
        current_chr = 1
        chr_start = 0
        offset = 0
        for number, line in enumerate(source):
            self.line_offsets.append(offset)
            offset += len(line)
            if not "chr" + str(current_chr) in line:
                self.chr_lines["chr" + str(current_chr)] = (chr_start, number)
                current_chr += 1
                chr_start = number

        print("Initialized")
        print(self.chr_lines)

    def readline(self, line):
        '''Quickly reads a given line from the MAF file by using the cached line line start position. The return value is a map from column name (as defined in MAF_Source.header) to value.'''
        self.source.seek(self.line_offsets[line])
        return {MAF_Source.header[i]: value for i, value in enumerate(self.source.readline().split())}

    def maf(self, chr, pos):
        '''Uses a binary search within the range of the given chromosome to quickly find the maf of the allele at the given position. If the allele does not exist, it returns -1'''
        lower, upper = self.chr_lines[chr]
        while lower < upper:
            midpoint = (upper - lower) // 2
            line = self.readline(midpoint)
            if line['position'] == pos: return line['maf']
            if line['position'] > pos: upper = midpoint
            else: lower = midpoint
        return -1

class Entry:
    '''A class that represents a single line entry of a snp stats file'''

    def __init__(self, data, maf_source=None):
        '''Parses a line of text into a snp stats line. Optionally, a MAF_Source object may be provided.'''
        self.data = data
        self.maf_source = maf_source

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
        return self.data[COLUMN['position']]

    def maf(self):
        '''Returns the MAF of the variants, either from the snp stats line or a separate file if provided.'''
        if self.maf_source: return self.maf_source.maf(self.chr(), self.pos())
        try:
            return float(self.data[COLUMN["minor_allele_frequency"]])
        except ValueError:
            return -1

    def structural_variant(self):
        '''Returns true if the variant is a structural variant. A structural variants is defined as any variant where either allele is not a single nucleotide'''
        return len(self.data[COLUMN["alleleA"]]) != 1 or len(self.data[COLUMN["alleleB"]]) != 1

    @staticmethod
    def for_entry(source, f, maf=None):
        '''Iterates over an entire snp stats file and runs f on each data line of the file'''
        entries = 0
        for line in source:
            if line[0] == "#": continue
            if HEADER in line: continue
            if ALT_HEADER in line: 
                use_alt_header()
                continue
            f(Entry(line.split(), maf))

    @staticmethod
    def for_folder(folder, f, maf=None):
        '''Iterates over an entire folder an runs f on each line in a snp stats file in the folder. The folder must have files for all 22 autosomes labeled chrN-snp-stats.txt'''
        for i in range(1, 23):
            Entry.for_entry(open(folder + '/' + file_name(i), 'r'), f, maf)


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
        self.n += 1
        self.average = self.average + (value - self.average)/self.n

    def in_range(self, item):
        return (self.minimum < item) and (item < self.maximum)

    def add_if_in_range(self, x, y):
        if self.in_range(x):
            self.add(y)
            return True
        else:
            return False

    def average(self):
        return self.average

    def count(self):
        return self.n

    @staticmethod
    def default_bins():
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
        for bin in bins:
            if bin.in_range(item): return bin

    @staticmethod
    def sort(bins, key, item):
        for bin in bins:
            if bin.add_if_in_range(key, item): return True
        return False

def write_header(f):
    f.write(HEADER + "\n")

def file_name(chrom):
    return "chr" + str(chrom) + "-snp-stats.txt"

