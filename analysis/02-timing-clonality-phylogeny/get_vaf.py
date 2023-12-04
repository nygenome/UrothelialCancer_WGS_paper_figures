#!/nfs/sw/python/python-3.7.1/bin/python
import pysam
import argparse
import pandas as pd
import itertools

# look for VAF at position in VCF

class Pileup():
    ''''''
    def __init__(self, pileupread, ref, alt,
                 variant_type):
        self.min_mq = 10
        self.min_bq = 10
        self.ref = ref
        self.alt = alt
        self.variant_type = variant_type
        self.pileupread = pileupread
        self.pos_in_read = pileupread.query_position
        self.read_name = pileupread.alignment.query_name
        self.is_ref = False
        self.is_alt = False
        self.get_count()
    
    def get_count(self):
        '''
        Variant calling for indels:

        For insertion:
            Check whether the length of the insertion and sequence matches alt allele.
            If there is no indel at the anchor position (even if the base at the anchor position doesn't match),
            we consider the read as adding support to the reference.
        For deletions:
            If the length of deletion matches alt allele, we consider that read supporting the
            alt allele.
            If there is no deletion at that position (even if there are mismatches in the bases spanning the deletion),
            it's considered to support reference allele.
        '''
        if self.check_span() and self.check_filters():
            if self.variant_type in ['SNV', 'MNV', 'COMPLEX']:
                read_ref = self.pileupread.alignment.query_sequence[self.pos_in_read:self.pos_in_read + len(self.ref)]
                read_alt = self.pileupread.alignment.query_sequence[self.pos_in_read:self.pos_in_read + len(self.alt)]
                if self.variant_type in ['SNV', 'MNV']:
                    if read_ref == self.ref:
                        self.is_ref = True
                    elif read_alt == self.alt:
                        self.is_alt = True
                elif self.variant_type in ['COMPLEX']:
                    if self.pileupread.indel == 0 and read_ref == self.ref:
                        self.is_ref = True
                    elif self.pileupread.indel == (len(self.alt) - len(self.ref)) and read_ref == self.alt:
                        self.is_alt = True
            elif self.variant_type in ['INS', 'DEL']:
                if self.pileupread.indel == 0:
                    self.is_ref = True
                elif self.variant_type in ['INS']:
                    read_alt = self.pileupread.alignment.query_sequence[self.pos_in_read+1:self.pos_in_read + len(self.alt)] == self.alt[1:]
                    if self.pileupread.indel == (len(self.alt) - len(self.ref)) and read_alt == self.alt:
                        self.is_alt = True
                elif self.variant_type in ['DEL']:
                    if self.pileupread.indel == (len(self.alt) - len(self.ref)):
                        self.is_alt = True

    def check_span(self):
        '''filter reads that don't span the indel'''
        if self.pos_in_read:
            if self.pos_in_read + len(self.ref) > self.pileupread.alignment.query_alignment_end \
                    or self.pos_in_read + len(self.alt) > self.pileupread.alignment.query_alignment_end:
                return False
        return True
        
    def check_filters(self):
        '''confirm read is above mq and bq and is a primary alignment'''
        # if the position in the read is .is_del pos is None so take next
        # skip reads where the position is already a deletion (is_del)
        if not self.pos_in_read:
            return False
        bq = self.pileupread.alignment.query_qualities[self.pos_in_read]
        mq = self.pileupread.alignment.mapping_quality
        if bq < self.min_bq \
                or mq < self.min_mq \
                or self.pileupread.alignment.is_supplementary:
            return False
        return True


class VCF_compare():
    
    def __init__(self, vcf_file, union_vcf,
                 t_bam_file, n_bam_file,
                 tumor, normal,
                 output_vcf,
                 max_indel_len_for_count=10):
        self.bcf_in = self.read_vcf(vcf_file)
        self.high_conf = self.read_union_vcf(union_vcf)
        self.bcf_out = self.write_vcf(output_vcf)
        self.t_bam = Bam(t_bam_file)
        self.n_bam = Bam(n_bam_file)
        self.tumor = tumor
        self.normal = normal
        self.max_indel_len_for_count = max_indel_len_for_count
        self.add_existing()
        self.discover_missing()
    
    def find_header(self, vcf_file):
        '''
            Get VCF header line numbers (because pandas can't skip
            based on > 1 character words. ## is a comment in VCF but
            # can occur in the VCF INFO fields.
            '''
        try:
            with open(vcf_file) as vcf:
                for i, line in enumerate(vcf):
                    if line.startswith('#'):
                        last = i
        except UnicodeDecodeError:
            with open(vcf_file, encoding='latin-1') as vcf:
                for i, line in enumerate(vcf):
                    if line.startswith('#'):
                        last = i
        return last + 1

    
    def load_utf_backup(self, vcf_file):
        '''
        Load searchable VCF for lines that are not ASCII.
        '''
        last = self.find_header(vcf_file)
        header = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'INFO']
        try:
            data = pd.read_csv(vcf_file, skiprows=last, sep='\t',
                               names=header, encoding = 'utf8', dtype = {'chrom' : str})
        except UnicodeDecodeError:
            data = pd.read_csv(vcf_file, skiprows=last, sep='\t',
                               names=header, encoding = 'latin-1', dtype = {'chrom' : str})
        return data
    
    def get_csqs(self, record):
        '''
            Get new INFO field results. Works with (python3). Skip lines starting with '#' for vcf backup.
            Because # can occur in info lines from VEP it can't be used as a comment character
        '''
        alt_count = len(record.alts)
        csq_values = []
        csq_dicts = {}
        for i in range(alt_count):
            try:
                csq_line = record.info[self.source][i]
            except UnicodeDecodeError: # for CSQ results with accents and other unexpected non-ascii characters (rare)
                line_data = vcf_backup[(self.high_conf_backup.chrom == record.chrom)
                                       & (self.high_conf_backup.pos == record.pos)
                                       & (self.high_conf_backup.ref == record.ref)
                                       & (self.high_conf_backup.alt == ','.join(record.alts))]
                csq_line = line_data.INFO.values.tolist()[0].split(self.source + '=')[1]
                csq_line = csq_line.split(';')[0]
            csq_line = csq_line.split(',')[i]
            csq_values = csq_line.split('|')
            try:
                # Works with (python2)
                csq_dict = dict(itertools.izip(self.csq_columns, csq_values))
            except AttributeError:
                # Works with (python3)
                csq_dict = dict(zip(self.csq_columns, csq_values))
            csq_dicts[i] = csq_dict
        return csq_dicts
    
    def write_vcf(self, output_vcf):
        bcf_out = pysam.VariantFile(output_vcf, 'w', header=self.bcf_in.header)
        return bcf_out
    
    def read_vcf(self, vcf_file):
        '''
            Read in annotated VCF file.
        '''
        bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
        return bcf_in
    
    def read_union_vcf(self, vcf_file):
        '''
            Read in first four columns of union "VCF" file.
        '''
        high_conf = pd.read_csv(vcf_file, sep='\t', comment='#',usecols=[0,1,3,4,7],
                                names=['CHROM', 'POS', 'REF', 'ALT', 'INFO'])
        return high_conf
    
    def is_too_long(self, ref, alt):
        '''
            Test if an INDEL or COMPLEX event is too long for computing allele counts using NYGC's pileup method given the length cut off.
            '''
        too_long = False
        if max([len(ref), len(alt)]) > self.max_indel_len_for_count:
            too_long = True
        return too_long
    
    def add_existing(self):
        '''write existing records'''
        for record in self.bcf_in.fetch():
            if record.info['HighConfidence']:
            #if not self.is_too_long(record.ref, record.alts[0]):
                if 'AD' in record.samples[self.normal].keys():
                    self.bcf_out.write(record)

    def discover_missing(self):
        '''add missing records'''
        for row in self.high_conf.iterrows():
            self.alt = row[1].ALT
            self.ref = row[1].REF
            self.pos = row[1].POS
            self.chrom = row[1].CHROM
            self.csq = row[1].INFO.split(';')[[idx for idx, s in enumerate(row[1].INFO.split(';')) if s.startswith('CSQ=')][0]].replace('CSQ=', '')
            # search for matching record in VCF
            if not self.is_too_long(self.ref, self.alt):
                for record in self.bcf_in.fetch(self.chrom, self.pos - 1, self.pos):
                    if not self.union_in_vcf(record) and not self.is_too_long(record.ref, record.alts[0]):
                        variant_type = record.info['TYPE']
                        self.t_bam.create_record(self.chrom, self.pos, self.ref, self.alt, variant_type)
                        self.n_bam.create_record(self.chrom, self.pos, self.ref, self.alt, variant_type)
                        # should this be pos + length?
                        if self.t_bam.vaf > 0:
                            try:
                                new_record = self.bcf_out.new_record(contig=self.chrom,
                                                                     alleles=(self.ref, self.alt),
                                                                     start=self.pos-1, stop=self.pos,
                                                                     filter='PASS')
                                new_record.samples[self.tumor]['DP'] = self.t_bam.dp
                                new_record.samples[self.normal]['DP'] = self.n_bam.dp
                                new_record.samples[self.tumor]['AD'] = (self.t_bam.ref_count, self.t_bam.alt_count)
                                new_record.samples[self.normal]['AD'] = (self.n_bam.ref_count, self.n_bam.alt_count)
                                new_record.info['TYPE'] = 'ADDED'
                                new_record.info['CSQ'] = self.csq
                                self.bcf_out.write(new_record)
                            except:
                                pass
#                            print(str(new_record))
#                        alleles=None, id=None, qual=None, filter=None, info=None, samples=None]

    def union_in_vcf(self, record):
        '''Check if a record is in the non-union VCF and high confidence'''
        if self.alt == record.alts[0] \
                and self.ref == record.ref \
                and record.info['HighConfidence']:
            return True
        return False


class Bam():
    
    def __init__(self, bam_file):
        self.min_mq = 10
        self.min_bq = 10
        self.bam_in = self.read_bam(bam_file)
    
    def read_bam(self, bam):
        bam_in = pysam.AlignmentFile(bam, "rb")
        return bam_in
    
    def create_record(self, chr, pos,  ref, alt, variant_type):
        '''Create a new record if the alts are found in the pileup'''
        refs, alts, others = self.pileup(self.bam_in, chr, pos, ref, alt, variant_type)
        self.ref_count, self.alt_count, self.dp = self.get_dp(refs, alts, others)
        if self.alt_count > 0:
            self.vaf = self.compute_vaf(self.alt_count, self.dp)
        else:
            self.vaf = 0

    def pileup(self, bam_in, chr, pos, ref, alt, variant_type):
        '''Run pileup for a variant'''
        refs = []
        alts = []
        others = []
        pileup = bam_in.pileup(chr, pos - 1, pos)
        for pileupcolumn in pileup:
            if pileupcolumn.pos == pos - 1:
                for pileupread in pileupcolumn.pileups:
                    read = Pileup(pileupread, ref, alt, variant_type)
                    if read.is_ref:
                        refs.append(read.read_name)
                    elif read.is_alt:
                        alts.append(read.read_name)
                    else:
                        others.append(read.read_name)
        return refs, alts, others

    def get_dp(self, refs, alts, others):
        '''Calc ref_count, alt_count, dp
            check sets to make sure reads don't show up in multiple sets
            supporting multiple calls
        '''
        set_ref_raw = set(refs)
        set_alt_raw = set(alts)
        set_other_raw = set(others)
        ref_reads_set = set_ref_raw - set_alt_raw - set_other_raw
        alt_reads_set = set_alt_raw - set_ref_raw - set_other_raw
        other_reads_set = set_other_raw - set_ref_raw - set_alt_raw
        all_reads_set = alt_reads_set|ref_reads_set|other_reads_set
        # tally set in ref and alt, non-ref/alt, all reads
        ref_count = len(ref_reads_set)
        alt_count = len(alt_reads_set)
        # other_count=len(other_reads_set) # not used
        dp = len(all_reads_set)
        return ref_count, alt_count, dp

    def compute_vaf(self, alt_count, dp):
        '''
        Compute VAF from dp and alt_count.
        '''
        return (0 if dp==0 else round(float(alt_count)/dp,4))



def main():
    parser = argparse.ArgumentParser(prog='add_nygc_allele_counts',
    description='Runs pileup on tumor and normal bam files to compute allele counts for bi-allelic SNV and Indel variants in VCF file and adds pileup format columns to the VCF file.', epilog='',
    formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=100, width=150))
    parser.add_argument('-t', '--tumor_bam', help = 'Tumor BAM file.', required=True)
    parser.add_argument('-n', '--normal_bam', help = 'Normal BAM file.', required=True)
    parser.add_argument('--tumor', help = 'Tumor.', required=True)
    parser.add_argument('--normal', help = 'Normal.', required=True)
    parser.add_argument('-v', '--vcf', help = 'VCF file.', required=True)
    parser.add_argument('-u', '--union-vcf',
                        help = ' '.join(['VCF file full of all filtered HighConfidence variants.',
                                         'Any variant in this file and missing in the VCF will be searched for in the BAM files.']),
                        required=True)
    parser.add_argument('-o', '--output', help = 'Output VCF file.', required=True)
    parser.add_argument('-i', '--max_indel_len_for_count',
                        help='Maximum indel or delin (complex event) length for generating counts',
                        default=10, type=int)
    args=parser.parse_args()

    VCF_compare(vcf_file=args.vcf,
                union_vcf= args.union_vcf,
                t_bam_file=args.tumor_bam,
                n_bam_file=args.normal_bam,
                tumor=args.tumor,
                normal=args.normal,
                output_vcf=args.output,
                max_indel_len_for_count=args.max_indel_len_for_count)

if __name__ == '__main__':
    main()
