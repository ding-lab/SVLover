// Process 1: QC and Mapping
process mapping {
    label 'mapping'
    input:
    path(normal_raw)
    path(tumor_raw)
	path(minimap_index)

    output:
    path("01_mapping/normal.sorted.bam")
    path("01_mapping/normal.sorted.bam.bai")
    path("01_mapping/tumor.sorted.bam")
	path("01_mapping/tumor.sorted.bam.bai")
    script:
    """
    mkdir -p 00_qc 01_mapping
    ${params.samtools} fastq -@ 10 -T'*' -0 normal.fq ${normal_raw}
	${params.minimap2} -y -Y -t 10 -ax map-hifi ${minimap_index} normal.fq >normal_tmp.sam 2> 01_mapping/normal.minimap2.log 
	${params.samtools} sort -@ 2 -o 01_mapping/normal.sorted.bam normal_tmp.sam
	rm normal_tmp.sam normal.fq
    ${params.samtools} index -@ 2 01_mapping/normal.sorted.bam
	${params.pandepth} -i 01_mapping/normal.sorted.bam -t ${params.cpu} -o 00_qc/normal.pandepth.tsv

	${params.samtools} fastq -@ 10 -T'*' -0 tumor.fq ${tumor_raw}
	${params.minimap2} -y -Y -t 10 -ax map-hifi ${minimap_index} tumor.fq  >tumor_tmp.sam 2> 01_mapping/tumor.minimap2.log 
	${params.samtools} sort -@ 2 -o 01_mapping/tumor.sorted.bam tumor_tmp.sam
	rm tumor_tmp.sam tumor.fq
    ${params.samtools} index -@ 2 01_mapping/tumor.sorted.bam
	${params.pandepth} -i 01_mapping/tumor.sorted.bam -t ${params.cpu} -o 00_qc/tumor.pandepth.tsv
    """
}


process germlineSV {
    label 'germlineSV'
    input:
    path(normalbam)
	path(repeat4pbsv)
	path(reference)

    output:
	path("02_GermlineSV/germlineSV.vcf")

    script:
    """
    mkdir -p 02_GermlineSV
    ${params.pbsv} discover ${normalbam} 02_GermlineSV/normal.svsig.gz --tandem-repeats ${repeat4pbsv} --hifi --sample normal
    ${params.pbsv} call -j 2 ${reference} 02_GermlineSV/normal.svsig.gz 02_GermlineSV/germlineSV.vcf

    """
}


process germlineSNV {
    label 'germlineSNV'
    input:
    path(normal_bam)
	path(normal_bamindex)
	path(reference)

    output:
	path("03_germlineSNV/deepvariant.vcf.gz")

    script:
    """
    mkdir 03_germlineSNV
	folder=\$(dirname "${params.reference}")
    run_deepvariant \\
        --model_type PACBIO \\
        --ref \${folder}/chm13v2.0.fa \\
        --reads ${normal_bam} \\
        --output_vcf 03_germlineSNV/deepvariant.vcf.gz \\
        --num_shards 5

    """
}


process phase {
    label 'phase'
    input:
    path(normalbam)
	path(tumorbam)
	path(normalbai)
	path(tumorbai)
	path(germlineSNV)
	path(germlineSV)
	path(reference)

    output:
	path("04_phase/normal_haplotagSV.bam")
	path("04_phase/normal_haplotagSV.bam.bai")
	path("04_phase/tumor_haplotagSV.bam")
	path("04_phase/tumor_haplotagSV.bam.bai")
	path("germlineSV.pass.phased.vcf.gz")
	path("phasesummary.tsv")

    script:
    """
    mkdir 04_phase
    ${params.bcftools} view -f 'PASS,.' ${germlineSNV} -O v -o deepvariant.pass.vcf
    sed -i -E '/^#CHROM/ s/([^\t]+)\$/default/' deepvariant.pass.vcf
	${params.bgzip} deepvariant.pass.vcf
	${params.tabix} deepvariant.pass.vcf.gz

	sed -i '/^#CHROM/s/normal/default/g' ${germlineSV}
	${params.bgzip} ${germlineSV}
	${params.tabix} ${germlineSV}.gz

	${params.bcftools} view -f 'PASS,.' ${germlineSV}.gz -O z -o germlineSV.pass.vcf.gz
	${params.tabix} germlineSV.pass.vcf.gz 

	${params.hiphase} --bam ${normalbam} --output-bam 04_phase/normal_haplotagSV.bam --bam ${tumorbam} --output-bam 04_phase/tumor_haplotagSV.bam --vcf deepvariant.pass.vcf.gz --output-vcf deepvariant.pass.phased.vcf.gz --vcf germlineSV.pass.vcf.gz  --output-vcf germlineSV.pass.phased.vcf.gz  --reference ${reference} --threads 5 --ignore-read-groups --blocks-file phaseblock.tsv --summary-file phasesummary.tsv
    """
}


process somaticSVnanomonsv {
    label 'somaticSVnanomonsv'

    input:
    path(normal_phase_bam)
	path(tumor_phase_bam)
	path(normal_phase_bai)
	path(tumor_phase_bai)
    path(reference)

    output:
    path("05_somaticSV/Nanomonsv/tumor.nanomonsv.result.corrected2.vcf")

    script:
    """
    folder=\$(dirname "${params.reference}")
    mkdir -p 05_somaticSV/Nanomonsv
    ${params.nanomonsv} parse ${tumor_phase_bam} 05_somaticSV/Nanomonsv/tumor
    ${params.nanomonsv} parse ${normal_phase_bam} 05_somaticSV/Nanomonsv/normal
    ${params.nanomonsv} get --control_prefix 05_somaticSV/Nanomonsv/normal --control_bam ${normal_phase_bam} --qv25 05_somaticSV/Nanomonsv/tumor ${tumor_phase_bam} \${folder}/chm13v2.0.fa
    
    awk '/^#/ { print; next } { if (/SVLEN/) { print; } else { gsub(/SVINSLEN/, \"SVLEN\"); print; } }'  05_somaticSV/Nanomonsv/tumor.nanomonsv.result.vcf |${params.bcftools} view -f 'PASS,.' | ${params.bcftools} sort  > 05_somaticSV/Nanomonsv/tumor.nanomonsv.result.corrected.vcf
    ${params.bcftools} annotate -x INFO/END 05_somaticSV/Nanomonsv/tumor.nanomonsv.result.corrected.vcf > 05_somaticSV/Nanomonsv/tumor.nanomonsv.result.corrected2.vcf
    rm 05_somaticSV/Nanomonsv/tumor.nanomonsv.result.corrected.vcf
    sed -i 's/SVTYPE=BND/SVTYPE=TRA/g' 05_somaticSV/Nanomonsv/tumor.nanomonsv.result.corrected2.vcf
    """
}

process somaticSVseverus {
    label 'somaticSVseverus'
    input:
    path(normal_phase_bam)
	path(tumor_phase_bam)
	path(normal_phase_bai)
	path(tumor_phase_bai)
	path(phasedSNV)
	path(repeat4pbsv)
    
    output:
    path("05_somaticSV/Severus/somatic_SVs/severus_somatic.vcf")

    script:
    """
    mkdir -p 05_somaticSV/Severus

    ${params.severus} --target-bam ${tumor_phase_bam} --control-bam ${normal_phase_bam} --out-dir 05_somaticSV/Severus -t 2  --phasing-vcf ${phasedSNV} --vntr-bed ${repeat4pbsv}  --use-supplementary-tag --output-read-ids
    if [ ! -e 05_somaticSV/Severus/somatic_SVs/severus_somatic.vcf ]; then
        echo "somatic_SVs VCF not found, check severus_all.vcf to see if there is no somatic SV"
        mkdir -p 05_somaticSV/Severus/somatic_SVs
        cp 05_somaticSV/Severus/all_SVs/severus_all.vcf 05_somaticSV/Severus/somatic_SVs/severus_somatic.vcf
    fi
	perl ${params.severuspolish} 05_somaticSV/Severus/somatic_SVs/severus_somatic.vcf 05_somaticSV/Severus/somatic_SVs/severus_somatic_corrected.vcf
	sed -i 's/SVTYPE=BND/SVTYPE=TRA/g' 05_somaticSV/Severus/somatic_SVs/severus_somatic_corrected.vcf

    """
}

process somaticSVsavana {
    label 'somaticSVsavana'
    input:
    path(normal_phase_bam)
	path(tumor_phase_bam)
	path(normal_phase_bai)
	path(tumor_phase_bai)

    output:
    path("05_somaticSV/SAVANA/somatic_savana.vcf")

    script:
    """
    mkdir -p 05_somaticSV/SAVANA
	folder=\$(dirname "${params.reference}")
    ${params.savana} run -t ${tumor_phase_bam} -n ${normal_phase_bam} --ref \${folder}/chm13v2.0.fa --ref_index \${folder}/chm13v2.0.fa.fai --outdir 05_somaticSV/SAVANA --length 50 --threads 1 --sample Tumor 
    ${params.savana} classify --vcf 05_somaticSV/SAVANA/Tumor.sv_breakpoints.vcf --pb --somatic_output 05_somaticSV/SAVANA/somatic_savana.vcf --output 05_somaticSV/SAVANA/classify.vcf 
    sed -i 's/SVTYPE=BND/SVTYPE=TRA/g' 05_somaticSV/SAVANA/somatic_savana.vcf
    """
}


process somaticSVsvisionpro {
    label 'somaticSVsvisionpro'
    input:
    path(normal_phase_bam)
	path(tumor_phase_bam)
	path(normal_phase_bai)
	path(tumor_phase_bai)
	path(reference)
    
    output:
    path("05_somaticSV/SVision2/somatic_svision.vcf")

    script:
    """
    mkdir -p 05_somaticSV/SVision2
	folder=\$(dirname "${params.reference}")
 	${params.svisionpro} --target_path ${tumor_phase_bam} --base_path ${normal_phase_bam} --genome_path \${folder}/chm13v2.0.fa --model_path ${params.svision_model} --out_path 05_somaticSV/SVision2 --sample_name Tumor --detect_mode somatic --preset hifi  --max_sv_size 500000 --skip_nearby
 	python ${params.extractop} --input_vcf 05_somaticSV/SVision2/*vcf  --extract somatic
 	mv 05_somaticSV/SVision2/*somatic*vcf 05_somaticSV/SVision2/somatic_tmp.vcf

 	awk 'BEGIN {OFS=\"\t\"}  /^#/ {print \$0} !/^#/ {if (\$8 ~ /SVTYPE=INS/) { gsub(/END=[^;]+;/, \"\", \$8)} print \$0}' 05_somaticSV/SVision2/somatic_tmp.vcf > 05_somaticSV/SVision2/somatic_svision.vcf
 	rm 05_somaticSV/SVision2/somatic_tmp.vcf
 	sed -i 's/SVTYPE=BND/SVTYPE=TRA/g' 05_somaticSV/SVision2/somatic_svision.vcf
    """
}


process mergeSomaticSV {
    label 'mergeSomaticSV'
    input:
    path(nanomonVCF)
	path(severusVCF)
	path(savanaVCF)
	path(svisionVCF)
    
    output:
    path("merged.2caller.vcf")

    script:
    """
    echo ${nanomonVCF} > somatic.vcflist
    echo ${svisionVCF} >> somatic.vcflist
    echo ${savanaVCF} >> somatic.vcflist
    echo ${severusVCF} >> somatic.vcflist
	${params.survivor} merge somatic.vcflist 1000 1 0 0 0 50 merged.2caller.vcf
	
    """
}










workflow {
    // Process 1: QC and Mapping
    normal_raw = Channel.fromPath(params.normalinput)
    tumor_raw = Channel.fromPath(params.tumorinput)
	minimap_index = Channel.fromPath(params.index)

    (normalbam,normalbai,tumorbam,tumorbai) = mapping(normal_raw, tumor_raw, minimap_index)

    // Process 2: germline SV
    repeat4pbsv = Channel.fromPath(params.repeat4pbsv)
    reference = Channel.fromPath(params.reference)
    germlineSV = germlineSV(normalbam, repeat4pbsv, reference)    
	
	// Process 3: germline SNV
	reference = Channel.fromPath(params.reference)
	germlineSNV = germlineSNV(normalbam, normalbai, reference)

	// Process 4: phase
	(normal_phase_bam, normal_phase_bai, tumor_phase_bam, tumor_phase_bai, phasedSNV, phasesummary) = phase(normalbam, tumorbam, normalbai, tumorbai, germlineSNV, germlineSV, reference)

	// Process 5: somatic SV
	nanomonVCF = somaticSVnanomonsv(normal_phase_bam, tumor_phase_bam, normal_phase_bai, tumor_phase_bai, reference)
	severusVCF = somaticSVseverus(normal_phase_bam, tumor_phase_bam, normal_phase_bai, tumor_phase_bai, phasedSNV, repeat4pbsv)
	savanaVCF = somaticSVsavana(normal_phase_bam, tumor_phase_bam, normal_phase_bai, tumor_phase_bai)
	svisionVCF = somaticSVsvisionpro(normal_phase_bam, tumor_phase_bam, normal_phase_bai, tumor_phase_bai, reference)
	SurvivorOutput = mergeSomaticSV(nanomonVCF, severusVCF, savanaVCF, svisionVCF)



}

