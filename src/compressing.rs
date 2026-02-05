use std::path::PathBuf;

use aho_corasick::AhoCorasick;
use anyhow::Result;

use deepsize::DeepSizeOf;

use crate::{
    helpers::{utils::NodeId, DeterministicHashMap, NodeRegistry, ReverseNodeRegistry},
    parser::{bufreader_from_compressed, parse_walk_seq_plainly, path_from_line, ByteLineReader},
};

pub fn compress_remaining_file(
    file: &PathBuf,
    mut rules: DeterministicHashMap<NodeId, Vec<NodeId>>,
    compressed_paths: &[bool],
    node_reg: &NodeRegistry,
    rev_reg: &ReverseNodeRegistry,
) -> Result<()> {
    unroll_rules(&mut rules);
    log::info!("Unrolled all rules");
    let mut lookup_table: DeterministicHashMap<Vec<NodeId>, NodeId> =
        rules.into_iter().map(|(k, v)| (v, k)).collect();
    lookup_table.shrink_to_fit();
    log::info!("Created lookup table of size: {}", lookup_table.deep_size_of() as f64 / 1_000_000_000.0);
    let ac = {
        let mut patterns: Vec<String> = Vec::with_capacity(2 * lookup_table.len());
        patterns.extend(lookup_table
            .keys()
            .map(|v| v.iter().map(|n| (*n).to_string()).collect())
        );
        // Add reverse-complement patterns
        patterns.extend(
            lookup_table
                .keys()
                .map(|v| v.iter().rev().map(|n| n.flip().to_string()).collect()),
        );
        let ac = AhoCorasick::builder()
            .match_kind(aho_corasick::MatchKind::LeftmostLongest)
            .build(&patterns)
            .expect("Can build Aho-Corasick structure");
        ac
    };
    log::info!("Done with setting up Aho-Corasick data structures of size {} bytes ({} Gb)", ac.memory_usage(), ac.memory_usage() as f64 / 1_000_000_000.0);

    let data = bufreader_from_compressed(file)?;
    let line_reader = ByteLineReader::new(data);
    let mut path_index = 0;
    line_reader.for_each(|line| {
        let h = path_from_line(&line, node_reg);
        if let Some((name, sequence)) = h {
            if !compressed_paths[path_index] {
                let compressed_text = compress_haplotype(&sequence, &lookup_table, &ac);
                let compressed_text: String = compressed_text
                    .into_iter()
                    .map(|n| rev_reg.get_directed_name(n))
                    .collect();
                println!("W\t{}\t{}", name, compressed_text);
            }
            path_index += 1;
        }
    });
    log::info!("Done with compressing remaining file");
    Ok(())
}

#[cfg(test)]
fn compress(
    haplotypes: &[Vec<NodeId>],
    mut rules: DeterministicHashMap<NodeId, Vec<NodeId>>,
) -> Vec<String> {
    unroll_rules(&mut rules);
    log::info!("Unrolled all rules");
    let mut solutions = Vec::new();
    let lookup_table: DeterministicHashMap<Vec<NodeId>, NodeId> =
        rules.into_iter().map(|(k, v)| (v, k)).collect();
    let mut patterns: Vec<String> = lookup_table
        .keys()
        .map(|v| v.iter().map(|n| (*n).to_string()).collect())
        .collect();
    // Add reverse-complement patterns
    patterns.extend(
        lookup_table
            .keys()
            .map(|v| v.iter().rev().map(|n| n.flip().to_string()).collect()),
    );
    let ac = AhoCorasick::builder()
        .match_kind(aho_corasick::MatchKind::LeftmostLongest)
        .build(&patterns)
        .expect("Can build Aho-Corasick structure");
    for haplotype in haplotypes {
        let text = compress_haplotype(haplotype, &lookup_table, &ac);
        let text = text.into_iter().map(|n| n.to_string()).collect();
        solutions.push(text);
    }
    solutions
}

// Unrolls rules so that no rule text contains meta-nodes
fn unroll_rules(
    rules: &mut DeterministicHashMap<NodeId, Vec<NodeId>>,
) {
    let keys: Vec<NodeId> = rules.keys().copied().collect();
    for meta_node in keys {
        unroll_rule(&meta_node, rules);
    }
}

fn unroll_rule(key: &NodeId, rules: &mut DeterministicHashMap<NodeId, Vec<NodeId>>) {
    if rules[key].iter().all(|n| !n.is_meta_node()) {
        return;
    }

    let mut new_seq: Vec<NodeId> = Vec::new();
    let old_seq = rules[key].clone();
    for node in old_seq {
        if node.is_meta_node() {
            unroll_rule(&node.get_forward(), rules);
            let seq = if node.is_forward() {
                rules[&node].clone()
            } else {
                rules[&node.get_forward()]
                    .iter()
                    .rev()
                    .map(|n| n.flip())
                    .collect()
            };
            new_seq.extend(seq.into_iter());
        } else {
            new_seq.push(node);
        }
    }
    *rules.get_mut(key).expect("Can find key meta-node") = new_seq;
}

fn compress_haplotype(
    haplotype: &[NodeId],
    lookup_table: &DeterministicHashMap<Vec<NodeId>, NodeId>,
    ac: &AhoCorasick,
) -> Vec<NodeId> {
    let haystack: String = haplotype.iter().map(|n| n.to_string()).collect();
    let mut matches = vec![];
    for mat in ac.find_iter(&haystack) {
        matches.push((mat.start(), mat.end()));
    }
    let mut text = haystack;
    for (start, end) in matches {
        let pattern_text = &text[start..end];
        let pattern_nodes = parse_walk_seq_plainly(pattern_text.as_bytes());
        let replacement = match lookup_table.get(&pattern_nodes) {
            Some(meta_node) => *meta_node,
            None => {
                // Lookup reverse-complement
                let rev_comp: Vec<NodeId> = pattern_nodes.iter().rev().map(|n| n.flip()).collect();
                lookup_table[&rev_comp].flip()
            }
        };
        // Fill replacement with spaces so replacement won't move
        // anything around
        if replacement.to_string().len() > end - start {
            panic!("Trying to fit in {} for {}", replacement, pattern_text);
        }
        let replacemen = format!("{: <1$}", replacement.to_string(), end - start);
        text.replace_range(start..end, &replacemen);
    }
    text.retain(|c| c != ' ');
    parse_walk_seq_plainly(text.as_bytes())
}

#[cfg(test)]
mod tests {
    use crate::helpers::utils::UndirectedNodeId;

    use super::*;

    fn get_node(id: u32, meta: bool, orientation: bool) -> NodeId {
        let undir = UndirectedNodeId(id, meta);
        NodeId(
            undir,
            if orientation {
                crate::helpers::utils::Orientation::Forward
            } else {
                crate::helpers::utils::Orientation::Backward
            },
        )
    }

    #[test]
    fn test_comress_simple() {
        let h = parse_walk_seq_plainly(">1>2>3>4>5".as_bytes());
        let mut rules = DeterministicHashMap::default();
        rules.insert(
            get_node(6, true, true),
            vec![get_node(1, false, true), get_node(2, false, true)],
        );
        rules.insert(
            get_node(7, true, true),
            vec![get_node(3, false, true), get_node(4, false, true)],
        );
        rules.insert(
            get_node(8, true, true),
            vec![
                get_node(6, true, true),
                get_node(7, true, true),
                get_node(5, false, true),
            ],
        );
        let res = compress(&vec![h], rules);
        assert_eq!(res, vec![">8"]);
    }
}
