workflow CompareVariants {
    File comparison_list
    #Map[String, File] gatk_vcfs
    #Map[String, File] gatk_index
    #Map[String, File] sv_vcfs
    File gatk_vcfs_file
    File gatk_index_file
    File sv_vcfs_file
    File easy_bed
    File medium_bed
    File hard_bed
    File ref_fasta

    Array[Array[String]] comparisons = read_tsv(comparison_list)
    Map gatk_vcfs = read_map(gatk_vcfs_file)
    Map gatk_index = read_map(gatk_index_file)
    Map sv_vcfs = read_map(sv_vcfs_file)

    scatter (comparison in comparisons) {

        String sample1=comparison[0]
        String sample2=comparison[1]
        call compare_GATK {
            input:
                ref_fasta=ref_fasta,
                truth_vcf=gatk_vcfs[sample1],
                vcf=gatk_vcfs[sample2],
                truth_vcf_index=gatk_index[sample1],
                vcf_index=gatk_index[sample2],
                easy_bed=easy_bed,
                medium_bed=medium_bed,
                hard_bed=hard_bed,
                output_prefix="${sample1}.${sample2}"
        }
        call compare_Lumpy {
            input:
                vcf1=sv_vcfs[sample1],
                vcf2=sv_vcfs[sample2],
                sample1=sample1,
                sample2=sample2
        }
    }
}

task compare_GATK {
    File ref_fasta
    File truth_vcf
    File vcf
    File truth_vcf_index
    File vcf_index
    File easy_bed
    File medium_bed
    File hard_bed
    String output_prefix
    command {
        echo "GiaB\t${easy_bed}" > stratification.tsv && \
        echo "difficult\t${hard_bed}" >> stratification.tsv && \
        echo "rest\t${medium_bed}" >> stratification.tsv && \
        HGREF=${ref_fasta} \
        /opt/hap.py/bin/hap.py \
        ${truth_vcf} \
        ${vcf} \
        -o ${output_prefix} \
        --stratification stratification.tsv \
        --location chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
        --threads 4 \
        --preprocess-truth
    }
    runtime {
        docker: pkrusche/hap.py
        cpu: "4"
        memory: "16 GB"
        disks: "local-disk " + 100 + " HDD"
        preemptible: 5
    }
    output {
        File output_vcf="${output_prefix}.vcf.gz"
        File output_summary="${output_prefix}.summary.csv"
        File output_extended="${output_prefix}.extended.csv"
    }
}

task compare_Lumpy {
    File vcf1
    File vcf2
    String sample1
    String sample2

    String bedpe1_name="${sample1}.bedpe"
    String bedpe2_name="${sample2}.bedpe"
    String output_name="${sample1}-${sample2}.counts.txt"
    String python="/opt/ccdg/python-2.7.12/bin/python"

    command {
        echo "Sample1	Sample2	Class	Count" > ${output_name} && \
        /opt/ccdg/python-2.7.12/bin/svtools vcftobedpe -i ${vcf1} -o ${bedpe1_name} && \
        /opt/ccdg/python-2.7.12/bin/svtools vcftobedpe -i ${vcf2} -o ${bedpe2_name} && \
        cat ${bedpe1_name} | perl -ape '$F[1] -= 1; $F[2]+=1; $F[4] -= 1; $F[5] += 1; $_ = join("\t", @F)."\n"' > ${bedpe1_name}.padded.bedpe && \
        cat ${bedpe2_name} | perl -ape '$F[1] -= 1; $F[2]+=1; $F[4] -= 1; $F[5] += 1; $_ = join("\t", @F)."\n"' > ${bedpe2_name}.padded.bedpe && \
        bedtools pairtopair -is -a ${bedpe1_name}.padded.bedpe -b ${bedpe2_name}.padded.bedpe -type both -slop 50 | sort -u | ${python} compare_single_sample_based_on_strand.py -l 0 | grep -v 'only' | sed "s/^/${sample1}	${sample2}	/" >> ${output_name} && \
        bedtools pairtopair -is -a ${bedpe1_name}.padded.bedpe -b ${bedpe2_name}.padded.bedpe -type notboth -slop 50 | sort -u | ${python} compare_single_sample_based_on_strand.py -l 0 | grep 'only' | sed "s/^/${sample1}	${sample2}	/" >> ${output_name} && \
        bedtools pairtopair -is -b ${bedpe1_name}.padded.bedpe -a ${bedpe2_name}.padded.bedpe -type notboth -slop 50 | sort -u | ${python} compare_single_sample_based_on_strand.py -l 0 | grep 'only' | sed "s/^/${sample1}	${sample2}	/" >> ${output_name}
    }
    runtime { #TODO need docker that includes bedtools and python script
        docker: "halllab/svtools@sha256:f2f3f9c788beb613bc26c858f897694cd6eaab450880c370bf0ef81d85bf8d45"
        preemptible: 5
    }
    output {
        File output_counts="${output_name}"
    }
}
