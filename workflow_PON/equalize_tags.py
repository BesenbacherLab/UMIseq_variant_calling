import pysam
import sys
from collections import defaultdict
import argparse
from os import remove

def read_pair_generator(bam):
    '''
    Generate read pairs from an indexed BAM file.
    Ignores reads not paired, is secondary or is supplementary
    '''

    read_pairs = defaultdict(lambda: [None, None])
    try:
        for read in bam.fetch(until_eof=True):
            if not read.is_paired or read.is_secondary or read.is_supplementary:
                continue
            read_name = read.query_name

            if read_name not in read_pairs:
                if read.is_read1:
                    read_pairs[read_name][0] = read
                else:
                    read_pairs[read_name][1] = read
            else:
                if read.is_read1:
                    yield read, read_pairs[read_name][1]
                else:
                    yield read_pairs[read_name][0], read
                del read_pairs[read_name]
         
    except OSError as e:
        print(e)


def main():
    description = f'''Equalize BAM tags that only exist on one read in pair.
      Resulting BAM file is written to STDOUT.
      Example: equalize_tags -q -t BX -t NI -T input.BAM > equalized.BAM'''
    parser = argparse.ArgumentParser(prog = "equalize_tags", description = description)
    parser.add_argument("-q", action = "store_true", default = False, help = "Add mate quality (MQ) tag")
    parser.add_argument("-t", "--tags", action = "append",  help = "Tag to equalize")
    parser.add_argument("-i", "--input", help = "Input BAM file")
    parser.add_argument("-o", "--output", help = "Output BAM file")
    o = parser.parse_args()

    input_bam = pysam.AlignmentFile(o.input, "rb", ignore_truncation = True)
    output_bam = pysam.AlignmentFile(o.output, "wb", template = input_bam)

    for read1, read2 in read_pair_generator(input_bam):
        if(o.q):
            read1.set_tag(tag = "MQ", value = read2.mapping_quality, value_type = "i", replace = False)
            read2.set_tag(tag = "MQ", value = read1.mapping_quality, value_type = "i", replace = False)

        for tag in o.tags:
            if read1.has_tag(tag):
                read2.set_tag(tag = tag, value = read1.get_tag(tag), replace = False)
            elif read2.has_tag(tag):
                read1.set_tag(tag = tag, value = read2.get_tag(tag), replace = False)

        output_bam.write(read1)
        output_bam.write(read2)
        

    input_bam.close()
    output_bam.close()


if __name__ == "__main__":
    main()

