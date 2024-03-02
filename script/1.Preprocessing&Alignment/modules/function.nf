
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.id
    meta.sample       = row.sample
    meta.sample_type  = row.sample_type
    meta.patient      = row.patient
    meta.platform     = row.platform
    meta.lane         = row.lane

    def array = []
    if (!file(row.forward_file_name).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.forward_file_name}"
    }
    if (!file(row.reverse_file_name).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.reverse_file_name}"
    }
    array = [ meta, [ file(row.forward_file_name), file(row.reverse_file_name) ] ]
    return array
}

def create_somatic_bam_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.id
    meta.sample       = row.sample
    meta.sample_type  = row.sample_type
    meta.patient      = row.patient
    meta.platform     = row.platform

    def array = []
    if (!file(row.tumor_bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> bam file does not exist!\n${row.tumor_bam}"
    }
    if (!file(row.tumor_bai).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> bai file does not exist!\n${row.tumor_bai}"
    }
    if (row.normal_bam) {
        if (!file(row.normal_bam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> bam file does not exist!\n${row.normal_bam}"
        }
        if (!file(row.normal_bai).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> bai file does not exist!\n${row.tumor_bai}"
        }
        array = [ meta, [file(row.tumor_bam), file(row.tumor_bai), file(row.normal_bam), file(row.normal_bai)] ]
    } else {
        array = [ meta, [file(row.tumor_bam), file(row.tumor_bai)] ]
    }

    return array
}

def create_OV_somatic_bam_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.tumor_sample
    meta.sample       = row.tumor_sample
    meta.normal       = row.normal_sample

    def array = []
    if (row.normal_bam) {
        array = [ meta, [file(row.tumor_bam), file(row.tumor_bai), file(row.normal_bam), file(row.normal_bai)] ]
    } else {
        array = [ meta, [file(row.tumor_bam), file(row.tumor_bai)] ]
    }

    return array
}


def create_maf_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample

    def array = []
    array = [ meta, file(row.maf) ]
    return array
}

def create_bam_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.id
    meta.sample       = row.sample
    meta.sample_type  = row.sample_type
    meta.patient      = row.patient
    meta.platform     = row.platform

    def array = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> bam file does not exist!\n${row.normal_bam}"
    }
    if (!file(row.bai).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> bai file does not exist!\n${row.tumor_bai}"
    }
    array = [ meta, file(row.bam), file(row.bai) ]
    return array
}

def create_vcf_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.tumor        = row.sample
    meta.normal       = row.normal

    def array = []
    array = [ meta, file(row.vcf), file(row.index) ]
    return array
}

def create_germ_vcf_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.sample           = row.sample

    def array = []
    array = [ meta, file(row.vcf), file(row.tbi) ]
    return array
}

def create_merge_vcf_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.tumor        = row.sample
    meta.normal       = row.normal

    def array = []
    array = [ meta, file(row.mutect_vcf), file(row.strelka_snv), file(row.strelka_indel), file(row.varscan_snv), file(row.varscan_indel), file(row.rescue_indel) ]
    return array
}


def create_noidx_vcf_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.tumor        = row.sample
    meta.normal       = row.normal

    def array = []
    array = [ meta, file(row.vcf) ]
    return array
}

def create_seg_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample

    def array = []
    array = [ meta, file(row.seg)]
    return array
}
def create_varscan_vcf_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.tumor        = row.sample
    meta.normal       = row.normal

    def array = []
    array = [ meta, file(row.snv), file(row.indel) ]
    return array
}

def create_absolute_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.tumor        = row.sample
    meta.normal       = row.normal

    def array = []
    array = [ meta, file(row.segement), file(row.maf) ]
    return array
}

def create_pyclone_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.id


    def array = []
    array = [ meta, file(row.c_input), file(row.ln_input) ]
    return array
}
