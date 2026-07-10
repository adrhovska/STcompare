#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/* parameters */

params.source_dir = null
params.reference_dir = null
params.sample_aligned = null
params.sample_reference = null

params.project_dir = "$PWD"
params.script_dir = "$baseDir"

params.py_env = "python_env"
params.r_env = "r_env"

params.counts1 = null
params.counts2 = null
params.spatial2 = null

/* Parameter checks */

if (!params.source_dir) {
    throw new IllegalArgumentException("Missing required parameter: --source_dir")
}
if (!params.reference_dir) {
    throw new IllegalArgumentException("Missing required parameter: --reference_dir")
}
if (!params.sample_aligned) {
    throw new IllegalArgumentException("Missing required parameter: --sample_aligned")
}
if (!params.sample_reference) {
    throw new IllegalArgumentException("Missing required parameter: --sample_reference")
}

/* output directories and folders */

align_pair_name = "${params.sample_aligned}_aligned_to_${params.sample_reference}"
landmark_pair_name = "${params.sample_aligned}_paired_to_${params.sample_reference}"

run_dir = "${params.project_dir}/STcompare_outputs/${align_pair_name}"

landmark_outdir = "${run_dir}/landmarks/${landmark_pair_name}"
stalign_outdir = "${run_dir}/STalign"
stcompare_outdir = "${run_dir}/STcompare"

/* input paths */

source_image = file("${params.source_dir}/spatial/tissue_hires_image.png", checkIfExists: true)
reference_image = file("${params.reference_dir}/spatial/tissue_hires_image.png", checkIfExists: true)

source_pos = file("${params.source_dir}/spatial/tissue_positions.csv", checkIfExists: true)
reference_pos = file("${params.reference_dir}/spatial/tissue_positions.csv", checkIfExists: true)

source_scale = file("${params.source_dir}/spatial/scalefactors_json.json", checkIfExists: true)
reference_scale = file("${params.reference_dir}/spatial/scalefactors_json.json", checkIfExists: true)

counts1 = file(params.counts1 ?: "${params.source_dir}/filtered_feature_bc_matrix.h5", checkIfExists: true)
counts2 = file(params.counts2 ?: "${params.reference_dir}/filtered_feature_bc_matrix.h5", checkIfExists: true)
spatial2 = file(params.spatial2 ?: "${params.reference_dir}/spatial", checkIfExists: true)


/* 1: LandmarkPicker.py */

process LANDMARK_PICKER {

    publishDir "${landmark_outdir}", mode: 'copy'

    input:
    path source_image
    path reference_image
    val sample_aligned
    val sample_reference

    output:
    path "out/${sample_aligned}_points.csv", emit: points1
    path "out/${sample_reference}_points.csv", emit: points2

    script:
    """
    source "\$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ${params.py_env}

    python "${params.script_dir}/LandmarkPicker.py" \\
      --image1 "${source_image}" \\
      --image2 "${reference_image}" \\
      --sample_aligned "${sample_aligned}" \\
      --sample_reference "${sample_reference}" \\
      --project_dir "."

    mkdir -p out

    cp "landmarks/${landmark_pair_name}/${sample_aligned}_points.csv" "out/${sample_aligned}_points.csv"
    cp "landmarks/${landmark_pair_name}/${sample_reference}_points.csv" "out/${sample_reference}_points.csv"
    """
}

/* 2: STalignCode.py
 * STalignCode.py imports qc_plots.py.
 * Therefore PYTHONPATH must include params.script_dir.
 */

process STALIGN {

    publishDir "${stalign_outdir}", mode: 'copy'

    input:
    path source_image
    path reference_image
    path source_pos
    path reference_pos
    path source_scale
    path reference_scale
    path points1
    path points2
    val sample_aligned
    val sample_reference

    output:
    path "STalign/*_barcodes.csv", emit: aligned_pos
    path "STalign/*", emit: stalign_outputs

    script:
    """
    source "\$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ${params.py_env}

    PYTHONPATH="${params.script_dir}:\${PYTHONPATH:-}" python "${params.script_dir}/STalignCode.py" \\
      --image1 "${source_image}" \\
      --image2 "${reference_image}" \\
      --pos1 "${source_pos}" \\
      --pos2 "${reference_pos}" \\
      --scale1 "${source_scale}" \\
      --scale2 "${reference_scale}" \\
      --points1 "${points1}" \\
      --points2 "${points2}" \\
      --sample_aligned "${sample_aligned}" \\
      --sample_reference "${sample_reference}" \\
      --project_dir "." \\
      --outdir "STalign" \\
      --alignment_method stalign
    """
}

/* 3: STcompare.R */

process STCOMPARE {

    publishDir "${stcompare_outdir}", mode: 'copy'

    input:
    path counts1
    path counts2
    path aligned_pos
    path spatial2
    val sample_aligned
    val sample_reference

    output:
    path "STcompare/*", emit: stcompare_outputs

    script:
    """
    source "\$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ${params.r_env}

    Rscript "${params.script_dir}/STcompare.R" \\
      --counts1 "${counts1}" \\
      --counts2 "${counts2}" \\
      --pos1 "${aligned_pos}" \\
      --spatial2 "${spatial2}" \\
      --outdir "STcompare" \\
      --sample_aligned "${sample_aligned}" \\
      --sample_reference "${sample_reference}"
    """
}

/* final workflow */

workflow {

    LANDMARK_PICKER(
        source_image,
        reference_image,
        params.sample_aligned,
        params.sample_reference
    )

    STALIGN(
        source_image,
        reference_image,
        source_pos,
        reference_pos,
        source_scale,
        reference_scale,
        LANDMARK_PICKER.out.points1,
        LANDMARK_PICKER.out.points2,
        params.sample_aligned,
        params.sample_reference
    )

    STCOMPARE(
        counts1,
        counts2,
        STALIGN.out.aligned_pos,
        spatial2,
        params.sample_aligned,
        params.sample_reference
    )

    emit:
    aligned_coordinates = STALIGN.out.aligned_pos
    stalign_outputs = STALIGN.out.stalign_outputs
    stcompare_outputs = STCOMPARE.out.stcompare_outputs
}