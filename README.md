# Long-read Tumor–Normal SV Pipeline (LOng read Variant ExploreR, SVLover)

A Nextflow pipeline for **tumor–normal long-read whole-genome sequencing (WGS)** analysis using **PacBio HiFi** data.  
The workflow performs mapping, germline variant calling, haplotype phasing, somatic structural variant (SV) detection with multiple callers, and SV merging.

---

## Workflow Overview

1. **Mapping & QC**  
   minimap2, samtools, pandepth

2. **Germline SV Calling**  
   pbsv

3. **Germline SNV Calling**  
   DeepVariant (PACBIO model)

4. **Phasing & Haplotagging**  
   HiPhase

5. **Somatic SV Calling**  
   - nanomonsv  
   - Severus  
   - SAVANA  
   - SVision-Pro  

6. **SV Merging**  
   SURVIVOR (SVs supported by ≥2 callers)

---

## Inputs

- Tumor HiFi BAM  
- Normal HiFi BAM  
- Reference genome (T2T-CHM13)  
- minimap2 index  
- Tandem repeat annotation (for pbsv / Severus)

---

## Outputs

- Phased tumor and normal BAM files  
- Germline SNV and SV VCFs  
- Somatic SV VCFs from individual callers  
- Merged somatic SV VCF (`merged.2caller.vcf`)

---

## Usage

```bash
nextflow run main.nf \
  --normalinput normal.bam \
  --tumorinput tumor.bam \
  --reference chm13v2.0.fa \
  --index minimap2_index
