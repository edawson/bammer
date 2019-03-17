task FreeBayesTask{
    File inputBAM
    File inputBAI
    File refFA
    File refFAI
    File inputVCFgz
    File inputVCFtbi
    File regionsBED

    Int diskGB

    Int? mincov = 15

    String outbase = basename(inputBAM, ".bam")
    String bedbase = basename(basename(regionsBED, ".gz"), ".bed")

    command{
        freebayes -@ ${inputVCFgz} \
        --report-all-haplotype-alleles \
        --targets ${regionsBED} \
        --report-monomorphic  \
        -f ${refFA} \
        --min-coverage ${mincov} \
        ${inputBAM} \
        > ${outbase}.${bedbase}.vcf
    }

    runtime{
        docker : "erictdawson/freebayes"
        cpu : 1
        memory : "6 GB"
        preemptible : 2
        disks : "local-disk " + diskGB + " HDD"
    }

    output{
        File FreeBayesVCF = "${outbase}.${bedbase}.vcf"
    }
}


task BedToRegionsTask{
    File bedFile
    Int diskGB

    String outbase = basename(basename(bedFile, ".gz"), ".bed")

    command{
        set -o pipefail
        cat ${bedFile} | bed_to_samtools_region > ${outbase}.regions.txt
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

task ShardBED{
    File inputBED

    String outbase = basename( basename(inputBED, ".gz"), ".bed")

    command{
        split -a 5 -d -l 10000 ${inputBED} ${outbase}.shard.bed.
    }

    runtime{
        docker : "erictdawson/base"
        cpu : 1
        memory : "1 GB" 
        disks : "local-disk 50 HDD"
    }

    output{
        Array[File] bedShards = glob("*.shard.bed*")
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
    File freebayesRegionsFile
    File refFA
    File refFAI
    File inputVCFgz
    File inputVCFtbi

    Int bigSZ = ceil(size(inputBAM, "GB") + size(refFA, "GB") + size(inputVCFgz, "GB")) + 100

    String outbase = basename(inputBAM, ".bam")

    Int brSZ = ceil(size(regionsBED, "GB") * 2) + 50

    call ShardBED{
        input:
            inputBED=regionsBED
    }

    
    scatter (reg in ShardBED.bedShards){
        call FreeBayesTask{
            input:
                inputBAM=inputBAM,
                inputBAI=inputBAI,
                refFA=refFA,
                refFAI=refFAI,
                inputVCFgz=inputVCFgz,
                inputVCFtbi=inputVCFtbi,
                regionsBED=reg,
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
