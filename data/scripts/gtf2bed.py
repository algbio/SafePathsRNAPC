# Changes a GFF/GTF file to a BED file by parsing the exons and creating BED format from them
# Author: Ariel Caceres (el.ariel.cl@gmail.com)

import optparse
import os
import sys


def parse_line(line):
    s = line.split('\t')
    return s[0], s[2], s[3], s[4], s[6], next(
        filter(lambda x: 'transcript_id' in x, s[8].split(';')), ' transcript_id ""'
    ).split(' ')[-1][1:-1]


def main(options):
    if not (options.input and options.exons and options.startCodons and options.stopCodons and options.cds and
            os.path.exists(options.input)):
        sys.stderr.write(f'Error in some of the parameters specified: {options}')
        sys.exit(1)

    gtf_file = open(options.input, 'r')
    transcripts = dict()
    for line in gtf_file:
        if line.startswith('#') or not line:
            continue

        chr_, type_, start, end, strand, transcript_id = parse_line(line)

        if type_ in ['exon', 'start_codon', 'stop_codon', 'CDS']:
            if transcript_id not in transcripts:
                transcripts[transcript_id] = {
                    'exons': list(),
                    'start_codons': list(),
                    'stop_codons': list(),
                    'CDSs': list()
                }
            transcripts[transcript_id][f'{type_}s'].append({
                'chr': chr_,
                'strand': strand,
                'start': int(start),
                'end': int(end)
            })

    gtf_file.close()

    exons_bed = open(options.exons, 'w')
    start_codons_bed = open(options.startCodons, 'w')
    stop_codons_bed = open(options.stopCodons, 'w')
    cds_bed = open(options.cds, 'w')

    for transcript_id in transcripts:

        exons = transcripts[transcript_id]['exons']
        start_codons = transcripts[transcript_id]['start_codons']
        stop_codons = transcripts[transcript_id]['stop_codons']
        CDSs = transcripts[transcript_id]['CDSs']

        if exons:
            first_exon = exons[0]
            last_exon = exons[-1]
            transcript_line = (
                f'{first_exon["chr"]}\t{first_exon["start"] - 1}\t{last_exon["end"]}\t{transcript_id}\t0\t'
                f'{first_exon["strand"]}\t{first_exon["start"] - 1}\t{last_exon["end"]}\t0\t{len(exons)}\t'
                f'{",".join(list(map(lambda exon: str(exon["end"] - exon["start"] + 1), exons)))}\t'
                f'{",".join(list(map(lambda exon: str(exon["start"] - first_exon["start"]), exons)))}\n'
            )

            exons_bed.write(transcript_line)

        if start_codons:
            first_start_codon = start_codons[0]
            last_start_codon = start_codons[-1]
            transcript_line = (
                f'{first_start_codon["chr"]}\t{first_start_codon["start"] - 1}\t{last_start_codon["end"]}\t'
                f'{transcript_id}\t0\t{first_start_codon["strand"]}\t{first_start_codon["start"] - 1}\t'
                f'{last_start_codon["end"]}\t0\t{len(start_codons)}\t'
                f'{",".join(list(map(lambda codon: str(codon["end"] - codon["start"] + 1), start_codons)))}\t'
                f'{",".join(list(map(lambda codon: str(codon["start"] - first_start_codon["start"]), start_codons)))}\n'
            )

            start_codons_bed.write(transcript_line)

        if stop_codons:
            first_stop_codon = stop_codons[0]
            last_stop_codon = stop_codons[-1]
            transcript_line = (
                f'{first_stop_codon["chr"]}\t{first_stop_codon["start"] - 1}\t{last_stop_codon["end"]}\t'
                f'{transcript_id}\t0\t{first_stop_codon["strand"]}\t{first_stop_codon["start"] - 1}\t'
                f'{last_stop_codon["end"]}\t0\t{len(stop_codons)}\t'
                f'{",".join(list(map(lambda codon: str(codon["end"] - codon["start"] + 1), stop_codons)))}\t'
                f'{",".join(list(map(lambda codon: str(codon["start"] - first_stop_codon["start"]), stop_codons)))}\n'
            )

            stop_codons_bed.write(transcript_line)

        if CDSs:
            first_CDS = CDSs[0]
            last_CDS = CDSs[-1]
            transcript_line = (
                f'{first_CDS["chr"]}\t{first_CDS["start"] - 1}\t{last_CDS["end"]}\t'
                f'{transcript_id}\t0\t{first_CDS["strand"]}\t{first_CDS["start"] - 1}\t'
                f'{last_CDS["end"]}\t0\t{len(CDSs)}\t'
                f'{",".join(list(map(lambda CDS: str(CDS["end"] - CDS["start"] + 1), CDSs)))}\t'
                f'{",".join(list(map(lambda CDS: str(CDS["start"] - first_CDS["start"]), CDSs)))}\n'
            )

            cds_bed.write(transcript_line)

    exons_bed.close()
    start_codons_bed.close()
    stop_codons_bed.close()
    cds_bed.close()


if __name__ == "__main__":
    parser = optparse.OptionParser(description="gtf2bed.py")
    parser.add_option("--input", "-i", help="input file (GFF/GTF)")
    parser.add_option("--exons", "-e", help="output (exons) file (BED)")
    parser.add_option("--startCodons", "-s", help="output (start codons) file (BED)")
    parser.add_option("--stopCodons", "-t", help="output (stop codons) file (BED)")
    parser.add_option("--cds", "-c", help="output (CDS) file (BED)")

    options, _ = parser.parse_args()

    main(options)
