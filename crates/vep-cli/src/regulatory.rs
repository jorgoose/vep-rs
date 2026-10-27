// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Cache-backed regulatory overlaps (Ensembl VEP 116.2).

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock};

use anyhow::Context;
use serde::de::{IgnoredAny, MapAccess, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer};
use vep_core::cache_manifest::CacheManifest;
use vep_core::consequence::{Consequence, FeatureType, TranscriptConsequence};
use vep_core::variant::{InputVariant, VariantClass};

use crate::cache_resolve::ResolvedCache;
use crate::config::Config;

#[derive(Deserialize)]
struct Feature {
    #[serde(deserialize_with = "crate::json_cache::deserialize_arc_str")]
    stable_id: Arc<str>,
    start: u64,
    end: u64,
    strand: i8,
    feature_type: FeatureType,
    #[serde(
        default,
        deserialize_with = "crate::json_cache::deserialize_optional_arc_str"
    )]
    biotype: Option<Arc<str>>,
    #[serde(default)]
    binding_matrix: Option<Matrix>,
}

#[derive(Deserialize)]
struct Matrix {
    stable_id: String,
    #[serde(default)]
    transcription_factors: Vec<String>,
    // Scores are not computed yet. Count positions without allocating the matrix.
    #[serde(default, rename = "frequencies", deserialize_with = "matrix_length")]
    length: u64,
}

fn matrix_length<'de, D: Deserializer<'de>>(d: D) -> Result<u64, D::Error> {
    struct Length;
    impl<'de> Visitor<'de> for Length {
        type Value = u64;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("a position-indexed binding matrix")
        }
        fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<u64, A::Error> {
            let mut length = 0;
            while map.next_entry::<IgnoredAny, IgnoredAny>()?.is_some() {
                length += 1;
            }
            Ok(length)
        }
        fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<u64, A::Error> {
            let mut length = 0;
            while seq.next_element::<IgnoredAny>()?.is_some() {
                length += 1;
            }
            Ok(length)
        }
    }
    d.deserialize_any(Length)
}

struct Chromosome {
    shards: Vec<PathBuf>,
    index: OnceLock<Index>,
}

struct Index {
    features: Vec<Feature>,
    // Maximum end through each start-sorted prefix; bounds backward overlap searches.
    max_end: Vec<u64>,
}

/// Regulatory JSON shards, enumerated eagerly and parsed once per queried chromosome.
pub struct RegulatorySource {
    chromosomes: HashMap<String, Chromosome>,
    manifest: Arc<CacheManifest>,
    loading: Mutex<()>,
}

impl RegulatorySource {
    pub(crate) fn open(
        config: &Config,
        cache: &ResolvedCache,
    ) -> anyhow::Result<Option<Arc<Self>>> {
        if !config.regulatory {
            return Ok(None);
        }
        let root = std::path::Path::new(&cache.dir).join("regulatory");
        if !cache.info.regulatory || !root.is_dir() {
            crate::warnings::global().warn("--regulatory has no effect: no regulatory data is available in the cache, so no TFBS_* or regulatory_region_* consequence is emitted");
            return Ok(None);
        }
        let mut chromosomes = HashMap::new();
        for entry in std::fs::read_dir(&root)? {
            let entry = entry?;
            if !entry.file_type()?.is_dir()
                || entry.file_name().as_encoded_bytes().starts_with(b".")
            {
                continue;
            }
            let mut shards = Vec::new();
            for shard in std::fs::read_dir(entry.path())? {
                let path = shard?.path();
                if path.extension().is_some_and(|ext| ext == "json") {
                    shards.push(path);
                }
            }
            shards.sort();
            chromosomes.insert(
                entry.file_name().to_string_lossy().into_owned(),
                Chromosome {
                    shards,
                    index: OnceLock::new(),
                },
            );
        }
        Ok(Some(Arc::new(Self {
            chromosomes,
            manifest: Arc::clone(&cache.manifest),
            loading: Mutex::new(()),
        })))
    }

    pub(crate) fn prewarm(&self, chr: &str) -> anyhow::Result<()> {
        let Some(chromosome) = self.chromosomes.get(chr) else {
            return Ok(());
        };
        if chromosome.index.get().is_some() {
            return Ok(());
        }
        let _guard = self.loading.lock().unwrap_or_else(|e| e.into_inner());
        if chromosome.index.get().is_some() {
            return Ok(());
        }
        let mut features = Vec::new();
        let mut seen = HashSet::new();
        for path in &chromosome.shards {
            let bytes =
                std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
            self.manifest.verify_bytes(path, &bytes)?;
            let shard: Vec<Feature> = serde_json::from_slice(&bytes)
                .with_context(|| format!("parsing {}", path.display()))?;
            for feature in shard {
                if seen.insert((feature.feature_type, Arc::clone(&feature.stable_id))) {
                    features.push(feature);
                }
            }
        }
        features.sort_by_key(|feature| feature.start);
        let mut end = 0;
        let max_end = features
            .iter()
            .map(|feature| {
                end = end.max(feature.end);
                end
            })
            .collect();
        let _ = chromosome.index.set(Index { features, max_end });
        Ok(())
    }

    pub(crate) fn annotate(&self, variant: &mut InputVariant) {
        if variant.no_allele || variant.variant_class == VariantClass::ReferenceBlock {
            return;
        }
        let Some(index) = self
            .chromosomes
            .get(&variant.chr)
            .and_then(|c| c.index.get())
        else {
            return;
        };
        // Do not normalize insertions: Perl overlap requires both flanks inside the feature.
        let start = variant.start;
        let end = variant.end;
        let initial_rows = variant.transcript_consequences.len();
        let mut i = index
            .features
            .partition_point(|feature| feature.start <= end);
        while i > 0 {
            i -= 1;
            if index.max_end[i] < start {
                break;
            }
            let feature = &index.features[i];
            if feature.end < start {
                continue;
            }
            if feature.feature_type == FeatureType::MotifFeature && feature.binding_matrix.is_none()
            {
                continue;
            }
            let mut row = TranscriptConsequence {
                transcript_id: Arc::clone(&feature.stable_id),
                feature_type: feature.feature_type,
                // StructuralVariationOverlapAllele_to_output_hash does not call the
                // regulatory/motif sequence-allele factories, so it has no BIOTYPE.
                biotype: if variant.is_structural {
                    None
                } else {
                    feature.biotype.clone()
                },
                strand: feature.strand,
                feature_start: feature.start,
                feature_end: feature.end,
                ..Default::default()
            };
            // These predicates share tier 2: Perl emits ablation/amplification AND the overlap term.
            if start <= feature.start && end >= feature.end {
                let deletion = if variant.is_structural {
                    matches!(
                        variant.variant_class,
                        VariantClass::StructuralDeletion | VariantClass::MobileElementDeletion
                    )
                } else {
                    variant.alt_allele().len() < variant.ref_allele.len()
                        || variant.alt_allele() == b"-"
                };
                let amplification = variant.is_structural
                    && matches!(
                        variant.variant_class,
                        VariantClass::Duplication
                            | VariantClass::TandemDuplication
                            | VariantClass::TandemRepeat
                    );
                if deletion {
                    row.consequences
                        .push(if feature.feature_type == FeatureType::MotifFeature {
                            Consequence::TfbsAblation
                        } else {
                            Consequence::RegulatoryRegionAblation
                        });
                } else if amplification {
                    row.consequences
                        .push(if feature.feature_type == FeatureType::MotifFeature {
                            Consequence::TfbsAmplification
                        } else {
                            Consequence::RegulatoryRegionAmplification
                        });
                }
            }
            row.consequences
                .push(if feature.feature_type == FeatureType::MotifFeature {
                    Consequence::TfBindingSiteVariant
                } else {
                    Consequence::RegulatoryRegionVariant
                });
            row.impact = row.consequences[0].impact();
            if let Some(matrix) = feature
                .binding_matrix
                .as_ref()
                .filter(|_| !variant.is_structural)
            {
                row.plugin_data
                    .insert("MOTIF_NAME".into(), matrix.stable_id.clone());
                row.plugin_data.insert(
                    "TRANSCRIPTION_FACTORS".into(),
                    matrix.transcription_factors.join(","),
                );
                let relative = i128::from(start) - i128::from(feature.start) + 1;
                let position = if feature.strand < 0 {
                    i128::from(matrix.length) - relative + 1
                } else {
                    relative
                };
                if position >= 1 && position <= i128::from(feature.end - feature.start + 1) {
                    row.plugin_data
                        .insert("MOTIF_POS".into(), position.to_string());
                }
            }
            variant.transcript_consequences.push(row);
        }
        if initial_rows == 0 && !variant.transcript_consequences.is_empty() {
            variant.intergenic_overlap = true;
        }
        if let Some(term) = variant
            .transcript_consequences
            .iter()
            .flat_map(|r| &r.consequences)
            .copied()
            .min_by_key(|c| c.rank())
        {
            variant.most_severe_consequence = Some(term);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source() -> (tempfile::TempDir, Arc<RegulatorySource>) {
        let dir = tempfile::tempdir().unwrap();
        let chr = dir.path().join("regulatory/21");
        std::fs::create_dir_all(&chr).unwrap();
        let fixture = include_bytes!("../tests/fixtures/regulatory/1-1000000.json");
        std::fs::write(chr.join("1-1000000.json"), fixture).unwrap();
        // Overlapping shards repeat features; the loader must emit each feature once.
        std::fs::write(chr.join("1000001-2000000.json"), fixture).unwrap();
        let cache = ResolvedCache {
            dir: dir.path().to_string_lossy().into_owned(),
            info: vep_core::cache_info::CacheInfo {
                regulatory: true,
                ..Default::default()
            },
            manifest: Arc::new(CacheManifest::load(dir.path()).unwrap()),
            assembly: Some("GRCh38".into()),
            fasta: None,
            synonyms: Arc::new(Default::default()),
        };
        let config = Config::from_argv(["vep", "--regulatory"]).unwrap();
        let source = RegulatorySource::open(&config, &cache).unwrap().unwrap();
        (dir, source)
    }

    fn variant(start: u64, end: u64, reference: &[u8], alternate: &[u8]) -> InputVariant {
        InputVariant::new(
            "21".into(),
            start,
            end,
            reference.to_vec(),
            alternate.to_vec(),
        )
    }

    #[test]
    fn lazy_chromosome_loading_deduplicates_and_reports_corruption() {
        let (dir, source) = source();
        assert!(source.chromosomes["21"].index.get().is_none());
        source.prewarm("missing").unwrap();
        assert!(source.chromosomes["21"].index.get().is_none());
        source.prewarm("21").unwrap();
        let mut v = variant(106, 106, b"C", b"G");
        source.annotate(&mut v);
        assert_eq!(v.transcript_consequences.len(), 2);
        assert!(v.intergenic_overlap);
        assert_eq!(
            v.most_severe_consequence,
            Some(Consequence::TfBindingSiteVariant)
        );
        let motif = v
            .transcript_consequences
            .iter()
            .find(|r| r.feature_type == FeatureType::MotifFeature)
            .unwrap();
        assert_eq!(motif.plugin_data["MOTIF_NAME"], "ENSPFM_FORWARD");
        assert_eq!(motif.plugin_data["MOTIF_POS"], "2");
        assert_eq!(motif.plugin_data["TRANSCRIPTION_FACTORS"], "TF_A,TF_B");
        // An unread chromosome must report malformed JSON rather than silently miss overlaps.
        let (_cold_dir, cold) = self::source();
        std::fs::write(cold.chromosomes["21"].shards[0].clone(), b"[").unwrap();
        assert!(cold
            .prewarm("21")
            .unwrap_err()
            .to_string()
            .contains("parsing"));
        drop(dir);
    }

    #[test]
    fn insertion_boundaries_reverse_positions_and_missing_matrices() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        for (start, end, expected) in [(105, 104, false), (106, 105, true), (109, 108, false)] {
            let mut v = variant(start, end, b"-", b"A");
            source.annotate(&mut v);
            assert_eq!(
                v.transcript_consequences
                    .iter()
                    .any(|r| r.transcript_id.as_ref() == "ENSM_FORWARD"),
                expected
            );
        }
        let mut reverse = variant(115, 116, b"AC", b"-");
        source.annotate(&mut reverse);
        let row = reverse
            .transcript_consequences
            .iter()
            .find(|r| r.transcript_id.as_ref() == "ENSM_REVERSE")
            .unwrap();
        assert_eq!(row.plugin_data["MOTIF_POS"], "4"); // variant start, not end
        let mut before = variant(114, 116, b"ACC", b"-");
        source.annotate(&mut before);
        let row = before
            .transcript_consequences
            .iter()
            .find(|r| r.transcript_id.as_ref() == "ENSM_REVERSE")
            .unwrap();
        assert!(!row.plugin_data.contains_key("MOTIF_POS"));
        let mut missing = variant(135, 135, b"A", b"G");
        source.annotate(&mut missing);
        assert!(missing.transcript_consequences.is_empty());
    }

    #[test]
    fn complete_deletion_and_duplication_keep_generic_overlap_term() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        for (class, term) in [
            (
                VariantClass::StructuralDeletion,
                Consequence::RegulatoryRegionAblation,
            ),
            (
                VariantClass::Duplication,
                Consequence::RegulatoryRegionAmplification,
            ),
        ] {
            let mut v = variant(100, 120, b"A", b"<DEL>");
            v.is_structural = true;
            v.variant_class = class;
            source.annotate(&mut v);
            let row = v
                .transcript_consequences
                .iter()
                .find(|r| r.feature_type == FeatureType::RegulatoryFeature)
                .unwrap();
            assert_eq!(
                &row.consequences[..],
                &[term, Consequence::RegulatoryRegionVariant]
            );
            assert!(row.biotype.is_none());
            for row in &v.transcript_consequences {
                assert!(
                    row.plugin_data.is_empty(),
                    "structural rows omit motif extras"
                );
            }
        }
    }
}
