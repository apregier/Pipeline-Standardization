workflow CallVariants {
    File id_list
    File exclude_regions
    File ref_fasta
    File ref_fasta_index
    File ref_fasta_dict

    Array[Array[File]] inputSamples = read_tsv(id_list)
    scatter (sample in inputSamples) {
        call GATK {
            input:
                sample_name=sample[0],
                cram_file=sample[1],
                cram_index=sample[2],
                ref_fasta=ref_fasta,
                ref_fasta_index=ref_fasta_index,
                ref_fasta_dict=ref_fasta_dict,
		output_vcf_filename="gatk.vcf.gz"
        }
        call ExtractReads {
            input:
                input_file=sample[1],
                input_index=sample[2],
        }
        call Lumpy {
            input:
                splitter_file=ExtractReads.splitter_file,
                discordants_file=ExtractReads.discordant_file,
                splitter_index=ExtractReads.splitter_index,
                discordants_index=ExtractReads.discordant_index,
                exclude_regions=exclude_regions,
                cram_file=sample[1]
        }
        call SVTyper {
            input:
                input_vcf=Lumpy.output_vcf,
                cram_file=sample[1],
		cram_index=sample[2],
                sample_name=sample[0]
        }
    }
    call GatherFiles as GatherSVs {
        input:
		pairs=SVTyper.outPair,
		mapName="SVTyperMap"
    }
    call GatherFiles as GatherGATK {
    	input:
		pairs=GATK.outPair,
		mapName="GATKMap"
    }
    call GatherFiles as GatherGATKIndex {
    	input:
		pairs=GATK.outPairIndex,
		mapName="GATKIndexMap"
    }
    output {
        File sv_vcfs = GatherSVs.mapFile
        File gatk_vcfs = GatherGATK.mapFile
        File gatk_index = GatherGATKIndex.mapFile
    }
}

task GatherFiles {
	Array[Pair[String, File]] pairs
	String mapName
	command {
        	cat ${write_lines(pairs)} > ${mapName}
        }
        runtime {
            docker: "broadinstitute/gatk3:3.5-0"
            memory: "4 GB"
            cpu: "1"
            preemptible: 5
        }
	output {
		File mapFile="${mapName}"
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
        docker: "broadinstitute/gatk3:3.5-0"
        memory: "16 GB"
        cpu: "2"
        disks: "local-disk " + 50 + " HDD"
    }
    output {
        File variants = "${output_vcf_filename}"
        File variants_index = "${output_vcf_filename}.tbi"
	Pair[String, File] outPair = ("${sample_name}", variants)
	Pair[String, File] outPairIndex = ("${sample_name}", variants_index)
    }
}

task ExtractReads {
    File input_file
    File input_index
    String splitter_file_name="splitters.bam"
    String discordant_file_name="discordants.bam"

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
        docker: "halllab/extract-sv-reads@sha256:192090f72afaeaaafa104d50890b2fc23935c8dc98988a9b5c80ddf4ec50f70c"
        cpu: "4"
        memory: "4 GB"
        disks: "local-disk " + 100 + " HDD"
        preemptible: 5
    }
    output {
        File splitter_file = "${splitter_file_name}"
        File discordant_file = "${discordant_file_name}"
        File splitter_index = "${splitter_file_name}.bai"
        File discordant_index = "${discordant_file_name}.bai"
    }
}

task Lumpy {
    File splitter_file
    File discordants_file
    File splitter_index
    File discordants_index
    File cram_file
    File exclude_regions
    String output_vcf_name="lumpy.vcf"

    command {
        lumpyexpress \
            -P \
            -T lumpy.temp \
            -o ${output_vcf_name} \
            -B ${cram_file} \
            -S ${splitter_file} \
            -D ${discordants_file} \
            -x ${exclude_regions} \
            -k \
            -v
    }
    runtime {
        docker: "halllab/lumpy@sha256:59ce7551307a54087e57d5cec89b17511d910d1fe9fa3651c12357f0594dcb07"
        cpu: "4"
        memory: "12 GB"
        disks: "local-disk " + 100 + " HDD"
        preemptible: 5
    }
    output {
        File output_vcf="${output_vcf_name}"
    }
}

task SVTyper {
    File input_vcf
    File cram_file
    File cram_index
    String sample_name

    String json_output_file_name="output.json"
    String vcf_output_file_name="output.vcf"
    command {
        cat ${input_vcf} | svtyper \
            -B ${cram_file} \
            -l ${json_output_file_name} \
            > ${vcf_output_file_name}
    }
    runtime {
        docker: "halllab/svtyper@sha256:21d757e77dfc52fddeab94acd66b09a561771a7803f9581b8cca3467ab7ff94a"
        cpu: "2"
        memory: "12 GB"
        disks: "local-disk " + 100 + " HDD"
        preemtible: 5
    }
    output {
        File output_vcf="${vcf_output_file_name}"
        File output_json="${json_output_file_name}"
	Pair[String, File] outPair=(sample_name, output_vcf)
    }
}
