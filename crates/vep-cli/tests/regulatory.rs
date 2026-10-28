// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use serde_json::Value;
use std::fs;
use std::path::Path;
mod common;

fn cache(dir: &Path) -> std::path::PathBuf {
    common::regulatory_cache(dir)
}

fn annotate(dir: &Path, cache: &Path, format: &str, regulatory: bool) -> String {
    annotate_with_args(dir, cache, format, regulatory, &[])
}

fn annotate_with_args(
    dir: &Path,
    cache: &Path,
    format: &str,
    regulatory: bool,
    extra_args: &[&str],
) -> String {
    let output = dir.join(format!("out-{format}-{regulatory}"));
    let mut command = common::isolated_vep(dir);
    command
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--no_stats",
            "--no_fasta",
            "--force_overwrite",
            "--format",
            "vcf",
            "--output_format",
            if format == "default" { "vep" } else { format },
        ])
        .arg("--input_file")
        .arg(dir.join("in.vcf"))
        .arg("--output_file")
        .arg(&output);
    if regulatory {
        command.arg("--regulatory");
    }
    command.args(extra_args);
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!String::from_utf8_lossy(&result.stderr).contains("no regulatory data"));
    fs::read_to_string(output).unwrap()
}

#[test]
fn regulatory_rows_and_extras_reach_every_text_format() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(dir.path().join("in.vcf"), "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tforward\tC\tG\t.\t.\t.\n21\t115\treverse\tA\tG\t.\t.\t.\n").unwrap();
    let json = annotate(dir.path(), &cache, "json", true);
    let records: Vec<Value> = json
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    let first = &records[0];
    assert_eq!(
        first["regulatory_feature_consequences"][0]["regulatory_feature_id"],
        "ENSR_FIXTURE"
    );
    assert_eq!(
        first["regulatory_feature_consequences"][0]["biotype"],
        "enhancer"
    );
    assert!(first["regulatory_feature_consequences"][0]
        .get("strand")
        .is_none());
    let motif = &first["motif_feature_consequences"][0];
    assert_eq!(motif["motif_feature_id"], "ENSM_FORWARD");
    assert_eq!(motif["motif_name"], "ENSPFM_FORWARD");
    assert_eq!(motif["motif_pos"], 2);
    assert_eq!(
        motif["transcription_factors"],
        serde_json::json!(["TF_A", "TF_B"])
    );
    assert!(motif.get("motif_score_change").is_none());
    assert!(motif.get("high_inf_pos").is_none());
    assert_eq!(first["most_severe_consequence"], "TF_binding_site_variant");
    assert_eq!(records[1]["motif_feature_consequences"][0]["motif_pos"], 4);
    assert_eq!(
        records[1]["motif_feature_consequences"][0]["transcription_factors"],
        serde_json::json!([])
    );
    for format in ["default", "tab", "vcf"] {
        let text = annotate(dir.path(), &cache, format, true);
        let rows = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n");
        assert!(
            rows.contains("regulatory_region_variant"),
            "{format}: {rows}"
        );
        assert!(rows.contains("TF_binding_site_variant"), "{format}: {rows}");
        assert!(rows.contains("enhancer"), "{format}: {rows}");
        assert!(rows.contains("ENSPFM_FORWARD"), "{format}: {rows}");
        assert!(rows.contains("TF_A"), "{format}: {rows}");
    }
    let disabled = annotate(dir.path(), &cache, "json", false);
    for line in disabled.lines() {
        let record: Value = serde_json::from_str(line).unwrap();
        assert!(record.get("regulatory_feature_consequences").is_none());
        assert!(record.get("motif_feature_consequences").is_none());
        assert_eq!(record["most_severe_consequence"], "intergenic_variant");
    }
    // Metadata is authoritative even when a regulatory tree exists.
    fs::write(cache.join("info.json"), r#"{"assembly":"GRCh38"}"#).unwrap();
    fs::write(cache.join("regulatory/21/1-1000000.json"), "[").unwrap();
    let mut command = common::isolated_vep(dir.path());
    let result = command
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--regulatory",
            "--no_fasta",
            "--no_stats",
            "--json",
            "--input_file",
        ])
        .arg(dir.path().join("in.vcf"))
        .args(["--output_file", "STDOUT"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stderr)
        .contains("no regulatory data is available in the cache"));
    assert!(!String::from_utf8_lossy(&result.stdout).contains("motif_feature_consequences"));
}

#[test]
fn regulatory_and_motif_rows_use_dash_for_absent_gene() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tforward\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    for format in ["default", "tab"] {
        let text = annotate(dir.path(), &cache, format, true);
        let features: Vec<_> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|columns| matches!(columns[5], "RegulatoryFeature" | "MotifFeature"))
            .map(|columns| (columns[4], columns[3]))
            .collect();
        assert_eq!(
            features,
            [("ENSR_FIXTURE", "-"), ("ENSM_FORWARD", "-")],
            "{format}: {text}"
        );
    }
}

fn assert_regulatory_only_rows(extra_args: &[&str], regulatory: bool, intergenic: bool) {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t110\tregulatory_only\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    for format in ["default", "tab"] {
        let text = annotate_with_args(dir.path(), &cache, format, true, extra_args);
        let rows: Vec<_> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect();
        let mut expected = Vec::new();
        if regulatory {
            expected.push((
                "ENSR_FIXTURE",
                "RegulatoryFeature",
                "regulatory_region_variant",
            ));
        }
        if intergenic {
            expected.push(("-", "-", "intergenic_variant"));
        }
        assert_eq!(
            rows.iter()
                .map(|row| (row[4], row[5], row[6]))
                .collect::<Vec<_>>(),
            expected,
            "{format} {extra_args:?}: {text}"
        );
        for row in &rows {
            assert_eq!(row[3], "-");
            if row[5] == "RegulatoryFeature" {
                if format == "default" {
                    assert!(row[13].split(';').any(|item| item == "BIOTYPE=enhancer"));
                } else {
                    let header = text
                        .lines()
                        .find(|line| line.starts_with("#Uploaded_variation"))
                        .unwrap();
                    let biotype = header
                        .split('\t')
                        .position(|column| column == "BIOTYPE")
                        .unwrap();
                    assert_eq!(row[biotype], "enhancer");
                }
            }
        }
    }
    let text = annotate_with_args(dir.path(), &cache, "json", true, extra_args);
    let records: Vec<Value> = text
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(records.len(), 1);
    let record = &records[0];
    assert!(record.get("transcript_consequences").is_none());
    assert!(record.get("motif_feature_consequences").is_none());
    assert_eq!(
        record["most_severe_consequence"],
        if regulatory {
            "regulatory_region_variant"
        } else {
            "?"
        }
    );
    if regulatory {
        assert_eq!(
            record["regulatory_feature_consequences"],
            serde_json::json!([{
                "biotype": "enhancer",
                "consequence_terms": ["regulatory_region_variant"],
                "impact": "MODIFIER",
                "regulatory_feature_id": "ENSR_FIXTURE",
                "variant_allele": "G"
            }])
        );
    } else {
        assert!(record.get("regulatory_feature_consequences").is_none());
    }
    if intergenic {
        assert_eq!(
            record["intergenic_consequences"],
            serde_json::json!([{
                "consequence_terms": ["intergenic_variant"],
                "impact": "MODIFIER",
                "variant_allele": "G"
            }])
        );
    } else {
        assert!(record.get("intergenic_consequences").is_none());
    }
    let text = annotate_with_args(dir.path(), &cache, "vcf", true, extra_args);
    let rows: Vec<_> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(rows.len(), 1);
    let info = rows[0].split('\t').nth(7).unwrap();
    if regulatory {
        let consequences: Vec<_> = info.strip_prefix("CSQ=").unwrap().split(',').collect();
        assert_eq!(consequences.len(), 1 + usize::from(intergenic));
        let regulatory_row: Vec<_> = consequences[0].split('|').collect();
        assert_eq!(
            &regulatory_row[..8],
            &[
                "G",
                "regulatory_region_variant",
                "MODIFIER",
                "",
                "",
                "RegulatoryFeature",
                "ENSR_FIXTURE",
                "enhancer"
            ]
        );
        if intergenic {
            let intergenic_row: Vec<_> = consequences[1].split('|').collect();
            assert_eq!(
                &intergenic_row[..8],
                &["G", "intergenic_variant", "MODIFIER", "", "", "", "", ""]
            );
        }
    } else {
        assert_eq!(info, ".");
    }
}

#[test]
fn regulatory_only_retains_intergenic_row() {
    assert_regulatory_only_rows(&[], true, true);
}

#[test]
fn regulatory_only_pick_retains_regulatory_row() {
    assert_regulatory_only_rows(&["--pick"], true, false);
}

#[test]
fn regulatory_only_coding_only_removes_all_rows() {
    assert_regulatory_only_rows(&["--coding_only"], false, false);
}

#[test]
fn regulatory_only_canonical_retains_intergenic_and_regulatory_rows() {
    assert_regulatory_only_rows(&["--canonical"], true, true);
}

#[test]
fn regulatory_only_biotype_retains_intergenic_and_regulatory_rows() {
    assert_regulatory_only_rows(&["--biotype"], true, true);
}

#[test]
fn regulatory_only_no_intergenic_removes_all_rows() {
    assert_regulatory_only_rows(&["--no_intergenic"], false, false);
}

/// A corrupt regulatory shard fails the embedding annotator's construction, as a
/// corrupt transcript shard does, instead of annotating every variant of the
/// chromosome without its regulatory rows and reporting the failure as an event.
#[test]
fn embedding_annotator_refuses_a_corrupt_regulatory_shard_at_construction() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    let config = vep_cli::Config::from_argv([
        "vep",
        "--json_cache",
        cache.to_str().unwrap(),
        "--regulatory",
        "--no_fasta",
    ])
    .unwrap();
    let error = match vep_cli::Annotator::from_config(&config) {
        Err(error) => format!("{error:#}"),
        Ok(_) => panic!("a cache with a corrupt regulatory shard built an annotator"),
    };
    assert!(
        error.contains("parsing") && error.contains("regulatory/22/1-1000000.json"),
        "{error}"
    );
}

#[test]
fn embedding_batches_load_regulatory_features_and_apply_summary_filters() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    // The embedding annotator loads every chromosome when it is built.
    fs::remove_dir_all(cache.join("regulatory/22")).unwrap();
    let config = vep_cli::Config::from_argv([
        "vep",
        "--json_cache",
        cache.to_str().unwrap(),
        "--regulatory",
        "--no_fasta",
    ])
    .unwrap();
    let annotator = vep_cli::Annotator::from_config(&config).unwrap();
    let mut variants = vec![
        vep_cli::InputVariant::new("21".into(), 106, 106, b"C".to_vec(), b"G".to_vec()),
        vep_cli::InputVariant::new("21".into(), 115, 115, b"A".to_vec(), b"G".to_vec()),
    ];
    annotator.annotate_batch(&mut variants);
    assert_eq!(
        variants[0].most_severe_consequence,
        Some(vep_cli::Consequence::TfBindingSiteVariant)
    );
    assert!(variants[0].intergenic_overlap);
    assert!(variants[1]
        .transcript_consequences
        .iter()
        .any(|r| r.plugin_data.get("MOTIF_POS").is_some_and(|p| p == "4")));
    let mut summary = variants.clone();
    vep_cli::apply_to_batch(
        &mut summary,
        &vep_cli::FilterConfig {
            summary: true,
            ..Default::default()
        },
    );
    assert_eq!(
        summary[0].summary_terms.as_deref(),
        Some(
            &[
                vep_cli::Consequence::TfBindingSiteVariant,
                vep_cli::Consequence::RegulatoryRegionVariant,
                vep_cli::Consequence::IntergenicVariant,
            ][..]
        )
    );
    let mut coding = variants.clone();
    vep_cli::apply_to_batch(
        &mut coding,
        &vep_cli::FilterConfig {
            row_filters: vep_cli::RowFilters {
                coding_only: true,
                ..Default::default()
            },
            ..Default::default()
        },
    );
    assert!(coding
        .iter()
        .all(|v| v.transcript_consequences.is_empty() && v.most_severe_consequence.is_none()));
    vep_cli::apply_to_batch(
        &mut variants,
        &vep_cli::FilterConfig {
            row_filters: vep_cli::RowFilters {
                no_intergenic: true,
                ..Default::default()
            },
            ..Default::default()
        },
    );
    assert!(variants
        .iter()
        .all(|v| v.transcript_consequences.is_empty() && v.most_severe_consequence.is_none()));
}
