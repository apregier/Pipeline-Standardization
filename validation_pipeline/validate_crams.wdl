workflow ValidateCrams {
    File comparison_list
    File id_list

    Array[Array[File]] inputSamples = read_tsv(id_list)
    Array[Array[File]] comparisons = read_tsv(comparison_list)

    scatter (sample in inputSamples) {
        call GATK {
            input:
                sample_name=sample[0],
                cram_file=sample[1],
                cram_index=sample[2],
                ref_fasta=ref_fasta,
                ref_fasta_index=ref_fasta_index,
                ref_fasta_dict=ref_fasta_dict,
        }
        call ExtractReads {
            input:
                
        }
        call Lumpy {
            input:
        }
        call SvTyper {
            input: 
        }
    }
    scatter (comparison in comparisons) {
        call compare_GATK {
            input:
                
        }
        call compare_Lumpy {
            input:
        }
    }
}

task GATK {
    String sample_name
    File cram_file
    File cram_index
    File ref_fasta
    File ref_fasta_index
    File ref_fasta_dict
    String output_vcf_filename
    command {
        java -Xmx10g -Xms10g -jar /usr/GenomeAnalysisTK.jar \
            -T HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${cram_file} \
            -o ${output_vcf_filename} \
            -rf BadCigar \
            --genotyping_mode DISCOVERY \
            --standard_min_confidence_threshold_for_calling 30 \
            --standard_min_confidence_threshold_for_emitting 0
    }
    runtime {
        docker: broadinstitute/gatk3:3.5-0
        memory: "16 GB"
        cpu: "2"
        disks: "local-disk " + 50 + " HDD"
    }
    output {
        File variants = "${output_vcf_filename}"
        File variants_index = "${output_vcf_filename}.tbi"
    }
}

task ExtractReads {
    command {
        extract-sv-reads --input-threads 4 \
            -e \
            -r \
            -i ${input_file} \
            -s ${splitter_file_name} \
            -d ${discordant_file_name} \
            && samtools index ${splitter_file_name} \
            && samtools index ${discordant_file_name}
    }
    runtime {
        docker: halllab/extract-sv-reads@sha256:192090f72afaeaaafa104d50890b2fc23935c8dc98988a9b5c80ddf4ec50f70c
        cpu: "4"
        memory: "4 GB"
        disks: "local-disk " + 100 + " HDD"
    }
    output {
        File splitter_file = ${splitter_file_name}
        File discordant_file = ${discordant_file_name}
        File splitter_index = ${splitter_index_name}
        File discordant_index = ${discordant_index_name}
    }
}

task Lumpy {
    String output_vcf_name

    command {
        lumpyexpress \
            -P \
            -T lumpy.temp \
            -o ${output_vcf_name} \
            -B ${cram_file} \
            -S ${splitters_file} \
            -D ${discordants_file} \
            -x ${exclude_regions} \
            -k \
            -v
    }
    runtime {
        docker: halllab/lumpy@sha256:59ce7551307a54087e57d5cec89b17511d910d1fe9fa3651c12357f0594dcb07
        cpu: "4"
        memory: "12 GB"
        disks: "local-disk " + 100 + " HDD"
    }
    output {
    }
}

task SVTyper {
    command {
        zcat ${input_vcf} | svtyper \
            -B ${cram_file} \
            -l ${json_output_file_name} \
            > ${vcf_output_file_name}
    }
    runtime {
        docker: halllab/svtyper@sha256:21d757e77dfc52fddeab94acd66b09a561771a7803f9581b8cca3467ab7ff94a
        cpu: "2"
        memory: "12 GB"
        disks: "local-disk " + 100 + " HDD"
    }
    output {
    }
}

task compare_GATK {
    File ref_fasta
    File truth_vcf
    File vcf
    File truth_vcf_index
    File vcf_index
    File stratification_file

    command {
        HGREF=${ref_fasta} \
        /opt/hap.py/bin/hap.py \
        ${truth_vcf} \
        ${vcf} \
        -o ${output_prefix} \
        --stratification ${stratification_file} \
        --location chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
        --threads 4 \
        --preprocess-truth
    }
    runtime {
        docker: pkrusche/hap.py
        cpu: "4"
        memory: "16 GB"
        disks: "local-disk " + 100 + " HDD"
    }
    output {
    }
}

task compare_Lumpy {
    command {
    }
    output {
    }
}
