#!/usr/bin/env python
from __future__ import print_function

import argparse
import numpy
import collections
import os
from operator import attrgetter


from Bio import SeqIO

import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class Asm:
    """
    This class contains a scaffoldsfile and statistics calculated on the basis
    of the fasta description headers found in the file.
    """
    headers = "Scaffolds_filename\tN50\t #contigs>=N50 \
        \tcoverage_contigs>N50\tminlen\tmaxlen\tavglen\tmedianlen\tstdlen \
        \tmincov\tmaxcov\tavgcov\tmediancov\tstdcov\t"

    def __init__(self, scaffoldsfile):
        """
        Takes a scaffoldsfile and calculates a range of stats based on it,
        mostly through methods found in the class

        Parameters
        ----------
        scaffoldsfile: string
            relative name of a SPAdes assembly

        Attributes (calculated)
        ----------
        scaffoldsfile: str
            relative path name of a SSPAdes assembly
        len_cov: namedTuple("LenCov", "len cov")
            contains the lengths and the coverages for the assembly
        n50: int
            the calculated N50 for the assembly
        stats_contigs_above_n50: namedtuple("Stats", "no_contigs avg_coverage")
            contains the number of contigs that together are longer than N50
            together with the average coverage for these contigs
        length_stats: (min, max, avg, median, std)
            five number summary for the contig lengths
        coverage_stats: (min, max, avg, median, std)
            five number summary for the contig coverages
        """
        self.scaffoldsfile = scaffoldsfile
        self.len_cov = self._length_coverage()
        self.n50 = self._calculate_n50()
        self.stats_contigs_above_n50 = self._stats_contigs_above_n50()
        self.length_stats = five_number_stats([x.len for x in self.len_cov])
        self.coverage_stats = five_number_stats([x.cov for x in self.len_cov])

    def _length_coverage(self):
        """
        Opens a scaffoldsfile and processes the header to access the length
        lengths and coverages for the assembly
        :return: list
            list containing named tuple with length and coverage
        """
        LenCov = collections.namedtuple("LenCov", "len cov")
        all_lencov = []

        with open(self.scaffoldsfile) as fh:
            print("processing " + self.scaffoldsfile)
            for record in SeqIO.parse(fh, "fasta"):
                fields = record.id.split("_")
                lencov = LenCov(len=int(fields[3]), cov=float(fields[5]))
                all_lencov.append(lencov)
        return all_lencov

    def _calculate_n50(self):
        """
        Calculates the N50 for the lengths found for this assembly
        :return: int
            the N50 value for the assembly
        """
        lengths = [x.len for x in self.len_cov]
        genomesize = sum(lengths)
        halfsize = genomesize/2

        n50 = "NA"
        lengths.sort(reverse=True)
        for counter in range(1, len(lengths)+1):
            this_n50 = sum(lengths[0:counter])
            if this_n50 > halfsize:
                n50 = lengths[counter-1]
                break
        return n50

    def _stats_contigs_above_n50(self):
        """
        Calculates the fewest number of contigs whose lengths together becomes
        more than the N50
        :return: namedtuple("Stats", "no_contigs avg_coverage")
            contains the number of contigs that together are longer than N50
            together with the average coverage for these contigs
        """
        len_cov_sorted = sorted(self.len_cov, key=attrgetter("len"), reverse=True)
        above = []
        i = 0
        while i < len(len_cov_sorted):
            if len_cov_sorted[i].len >= self.n50:
                above.append(len_cov_sorted[i])
            i += 1

        Stats = collections.namedtuple("Stats", "no_contigs avg_coverage")
        no_contigs_above = len(above)
        avg_coverage_above = sum([x.cov for x in above])/len(above)
        stats = Stats(no_contigs=no_contigs_above, avg_coverage=avg_coverage_above)
        return stats



    def __str__(self):
        """
        Creates nicely formatted output for writing to file/screen
        :return: string
            well formatted output
        """
        outstring = self.scaffoldsfile
        outstring += "\t" + str(self.n50)
        outstring += "\t" + str(self.stats_contigs_above_n50.no_contigs)
        outstring += "\t" + str(self.stats_contigs_above_n50.avg_coverage)
        outstring += "\t" + "\t".join(map(str, self.length_stats))
        outstring += "\t" + "\t".join(map(str, self.coverage_stats))
        return outstring


def five_number_stats(numberlist):
    """
    Calculates five number statistics on a list of numbers
    :param numberlist:
        list containing numbers
    :return: (min, max, avg, median, std) - tuple with floats
        tuple containing five number statics
    """
    minval = numpy.min(numberlist)
    maxval = numpy.max(numberlist)
    avg = numpy.mean(numberlist)
    median = numpy.median(numberlist)
    std = numpy.std(numberlist)
    return (minval, maxval, avg, median, std)


def calculate_directory_stats(asm_list):
    """
    Takes a list of Asms and creates nicely formatted results for printing
    :param asm_list: string
        list of Asms
    :return: string
        output that can be printed to file
    """
    output = ""
    for asm in asm_list:
        output += "\n" + asm.__str__()
    return output


def process_directory(directory):
    """
    Examines a directory for SPAdes scaffolds.fasta files and creates Asm
    objects from them (thus calculating asm stats), and then getting the
    results for printing.
    :param directory: string
        directory containing assemblies
    :return: string
        text string containing output that will be written to file
    """
    asm_list = []
    for path, dirs, files in os.walk(directory):
        for filename in files:
            if filename == "scaffolds.fasta":
                fullname = os.path.join(path, filename)
                new_asm = Asm(fullname)
                asm_list.append(new_asm)
    print("Have collected %s number of asms" % len(asm_list))
    output = calculate_directory_stats(asm_list)
    return output


def print_output(results, prefix):
    """
    Pretty writes results to the output file
    :param results: string
        results to be written to file
    :param prefix: string
        string that will be prefixed to output files
    """
    fo = open(prefix + ".out", "w")
    fo.write(Asm.headers)
    fo.write(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", metavar="DIRECTORY",
                        help="Directory containing directories with " +
                             "scaffolds.fasta in them")
    parser.add_argument("-p", "--prefix", metavar="STRING",
                        help="Prefix string for output files")
    args = parser.parse_args()

    results = process_directory(args.directory)
    print_output(results, args.prefix)
