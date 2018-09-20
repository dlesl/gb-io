extern crate string_cache_codegen;

use std::env;
use std::path::Path;
use std::fs::File;

fn main() {
    let qualifiers = [
        "allele",
        "anticodon",
        "bound_moiety",
        "chromosome",
        "codon_recognized",
        "codon_start",
        "country",
        "db_xref",
        "EC_number",
        "estimated_length",
        "exception",
        "experiment",
        "function",
        "gap_type",
        "gene",
        "gene_synonym",
        "inference",
        "isolation_source",
        "label",
        "linkage_evidence",
        "locus_tag",
        "map",
        "mobile_element_type",
        "mol_type",
        "ncRNA_class",
        "nomenclature",
        "note",
        "number",
        "old_locus_tag",
        "organelle",
        "organism",
        "product",
        "protein_id",
        "pseudo",
        "recombination_class",
        "regulatory_class",
        "ribosomal_slippage",
        "rpt_type",
        "standard_name",
        "strain",
        "sub_strain",
        "tissue_type",
        "transcript_id",
        "translation",
        "transl_except",
        "transl_table",
    ];
    let feature_kinds = [
        "assembly_gap",
        "CDS",
        "centromere",
        "C_region",
        "D-loop",
        "D_segment",
        "exon",
        "gene",
        "J_segment",
        "misc_feature",
        "misc_recomb",
        "misc_RNA",
        "mobile_element",
        "mRNA",
        "ncRNA",
        "precursor_RNA",
        "protein_bind",
        "regulatory",
        "repeat_region",
        "rep_origin",
        "rRNA",
        "source",
        "tmRNA",
        "tRNA",
        "V_segment",
    ];
    let file = File::create(&Path::new(&env::var("OUT_DIR").unwrap()).join("atoms.rs")).unwrap();
    string_cache_codegen::AtomType::new("QualifierKey", "qualifier_key!")
        .atoms(qualifiers.iter())
        .write_to(&file)
        .unwrap();
    string_cache_codegen::AtomType::new("FeatureKind", "feature_kind!")
        .atoms(feature_kinds.iter())
        .write_to(&file)
        .unwrap();
}
