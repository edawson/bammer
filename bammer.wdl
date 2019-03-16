task FreeBayesTask{
    File inputBAM
    File inputBAI
    File refFA
    File refFAI
    File inputVCFgz
    File inputVCFtbi
    String region

    Int diskGB

    Int? mincov = 15

    String outbase = basename(inputBAM, ".bam")

    command{
        freebayes -@ ${inputVCFgz} \
        --report-all-haplotype-alleles \
        --report-monomorphic  \
        -f ${refFA} \
        --min-coverage ${mincov} \
        ${inputBAM} ${region} \
        > ${outbase}.${region}.vcf
    }

    runtime{
        docker : "erictdawson/freebayes"
        cpu : 1
        memory : "6 GB"
        preemptible : 2
        disks : "local-disk " + diskGB + " HDD"
    }

    output{
        File FreeBayesVCF = "${outbase}.${region}.vcf"
    }
}

task SplitBedByCoverage{
    File inputBAM
    File inputBAI
    File inputBED
    Int diskGB 

    String outbase = basename(basename(inputBED, ".gz"), ".bed")

    command{
        split_bed_by_index -m ${inputBED} ${inputBAM} ${outbase}.covRestricted.bed
    }

    runtime{
        docker : "erictdawson/bedtools"
        cpu : 1
        memory : "3 GB"
        preemptible : 2
        disks : "local-disk " + diskGB + " HDD"
    }

    output{
        File coverageSplitBed = "${outbase}.covRestricted.bed"
    }
}



task BedToRegionsTask{
    File bedFile
    Int diskGB

    String outbase = basename(basename(bedFile, ".gz"), ".bed")

    command{
        zcat ${bedFile} | bed_to_samtools_region > ${outbase}.regions.txt
    }
    runtime{
        docker : "erictdawson/bedtools"
        cpu : 1
        memory : "1 GB"
        preemptible : 2
        disks : "local-disk " + diskGB + " HDD"
    }

    output{
        File regionsFile = "${outbase}.regions.txt"
    }
}

task FreeBayesMergeTask{
    Array[File] inputVCFs
    String outbase
    Int diskGB

    command{
        cat_file_of_files ${write_lines(inputVCFs)} | \
        vcffirstheader | \
        vcfstreamsort -w 1000 | \
        vcfuniq > ${outbase}.vcf
    }

    runtime{
        docker : "erictdawson/freebayes"
        cpu : 2
        memory : "7 GB"
        preemptible : 2
        disks : "local-disk " + diskGB + " HDD"
    }

    output{
        File mergedVCF = "${outbase}.vcf"
    }
}

task BGZTBI{
    File inputVCF
    Int diskGB

    String outbase = basename( basename(inputVCF, ".gz"), ".vcf")

    command {
        bgzip -c ${inputVCF} > ${outbase}.vcf.gz && \
        tabix ${outbase}.vcf.gz
    }

    runtime {
        docker : "erictdawson/base"
        cpu : 1
        memory : "1.6 GB"
        preemptible : 2
        disks : "local-disk " + diskGB + " HDD" 
    }

    output{
        File vcfGZ = "${outbase}.vcf.gz"
        File vcfTBI = "${outbase}.vcf.gz.tbi"
    }
}

workflow FreeBayesForceCall{
    File inputBAM
    File inputBAI
    File regionsBED
    File refFA
    File refFAI
    File inputVCFgz
    File inputVCFtbi

    Int bigSZ = ceil(size(inputBAM, "GB") + size(refFA, "GB") + size(inputVCFgz, "GB")) + 100

    String outbase = basename(inputBAM, ".bam")

    call SplitBedByCoverage{
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            inputBED=regionsBED,
            diskGB=bigSZ
    }

    Int brSZ = ceil(size(SplitBedByCoverage.coverageSplitBed, "GB") * 2) + 50

    call BedToRegionsTask{
        input:
            bedFile=SplitBedByCoverage.coverageSplitBed,
            diskGB=brSZ
    }
    

    
    scatter (reg in read_lines(BedToRegionsTask.regionsFile)){
        call FreeBayesTask{
            input:
                inputBAM=inputBAM,
                inputBAI=inputBAI,
                refFA=refFA,
                refFAI=refFAI,
                inputVCFgz=inputVCFgz,
                inputVCFtbi=inputVCFtbi,
                region=reg,
                diskGB=bigSZ
        }
    }

    Int frSZ = 500

    call FreeBayesMergeTask{
        input:
            inputVCFs=FreeBayesTask.FreeBayesVCF,
            outbase=outbase,
            diskGB=frSZ
    }

    call BGZTBI{
        input:
            inputVCF=FreeBayesMergeTask.mergedVCF,
            diskGB=frSZ
    }
}
