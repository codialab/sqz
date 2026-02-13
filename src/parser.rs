use anyhow::{anyhow, Result};
use flate2::read::MultiGzDecoder;
use itertools::Itertools;
use lazy_static::lazy_static;
use rand::{Rng, SeedableRng};
use regex::bytes::Regex;
use std::{io::BufReader, str::FromStr};

use crate::{
    decoding::decode_walk,
    helpers::{
        utils::{Address, Digram, LocalizedDigram, NodeId, Orientation, UndirectedNodeId},
        DeterministicHashMap, DeterministicHashSet, NodeRegistry, PathSegment, ReverseNodeRegistry,
    },
};
use std::{
    io::{self, BufRead, Read},
    path::PathBuf,
};

lazy_static! {
    static ref RE_WALK: Regex = Regex::new(r"([><])([!-;=?-~]+)").unwrap();
}

type NamedPath = Vec<(PathSegment, Vec<LocalizedDigram>)>;

pub struct ByteLineReader<R: io::Read + Send> {
    data: io::BufReader<R>,
}

impl<R: io::Read + Send> ByteLineReader<R> {
    pub fn new(data: R) -> Self {
        Self {
            data: BufReader::new(data),
        }
    }
}

impl<R: io::Read + Send> Iterator for ByteLineReader<R> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = Vec::new();
        match self.data.read_until(b'\n', &mut buf) {
            Err(e) => {
                log::error!("Error while reading graph: {e:?} ");
                Some(buf)
            }
            Ok(1..) => Some(buf),
            Ok(0) => None,
        }
    }
}

pub fn bufreader_from_compressed(file: &PathBuf) -> Result<io::BufReader<Box<dyn Read + Send>>> {
    let f = std::fs::File::open(file)?;
    if let Some(name) = file.file_name() {
        let reader: Box<dyn Read + Send> = if name.to_str().unwrap().ends_with(".gz") {
            log::info!("assuming that {} is gzip compressed..", file.display());
            Box::new(MultiGzDecoder::new(f))
        } else {
            Box::new(f)
        };
        Ok(io::BufReader::new(reader))
    } else {
        Err(anyhow!("Filename {:?} is not a valid filename", file))
    }
}

pub fn compare_file(
    file: &PathBuf,
    walks: &[Vec<NodeId>],
    grammar: &Grammar,
    node_ids_by_name: &NodeRegistry,
) {
    let data = bufreader_from_compressed(file).expect("Can read file");
    compare_file_content(data, walks, grammar, node_ids_by_name);
}

fn compare_file_content<R: Read + Send>(
    reader: R,
    walks: &[Vec<NodeId>],
    grammar: &Grammar,
    node_ids_by_name: &NodeRegistry,
) {
    let line_reader = ByteLineReader::new(reader);
    let rev_node: ReverseNodeRegistry = node_ids_by_name.clone().into();
    let mut counter = 0;
    let mut found_error = false;
    line_reader.for_each(|line| {
        let h = path_from_line(&line, node_ids_by_name);
        if let Some((haplotype_name, expected)) = h {
            let actual = decode_walk(walks[counter].clone(), grammar);
            for i in 0..expected.len() {
                if i > actual.len() - 1 {
                    log::error!(
                        "Difference at {} ({}: {}): - (calculated) vs {} (expected)",
                        haplotype_name.to_path_string(),
                        counter,
                        i,
                        rev_node.get_directed_name(expected[i])
                    );
                } else if expected[i] != actual[i] {
                    log::error!(
                        "Difference at {} ({}: {}): {} (calculated) vs {} (expected)",
                        haplotype_name.to_path_string(),
                        counter,
                        i,
                        rev_node.get_directed_name(actual[i]),
                        rev_node.get_directed_name(expected[i])
                    );
                    found_error = true;
                }
            }
            counter += 1;
        }
    });

    if found_error {
        return;
    }

    println!("No differences were found!");
}

pub fn parse_file_to_haplotypes_with_grammar(
    file: &PathBuf,
    should_print_other_lines: bool,
) -> Result<(Vec<Path>, NodeRegistry, Grammar)> {
    let data = bufreader_from_compressed(file)?;
    let (node_ids_by_name, _) = parse_node_ids(data, true)?;
    let data = bufreader_from_compressed(file)?;
    let grammar = parse_grammar(data, &node_ids_by_name)?;
    let data = bufreader_from_compressed(file)?;
    let (paths, nodes) =
        parse_file_content_to_haplotypes(data, node_ids_by_name, should_print_other_lines)?;
    Ok((paths, nodes, grammar))
}

pub fn parse_file_to_digrams(
    file: &PathBuf,
    should_print_other_lines: bool,
) -> Result<(NamedPath, NodeRegistry, DeterministicHashMap<usize, NodeId>)> {
    let data = bufreader_from_compressed(file)?;
    let (node_ids_by_name, _) = parse_node_ids(data, false)?;
    let data = bufreader_from_compressed(file)?;
    parse_file_content_to_digrams(data, node_ids_by_name, should_print_other_lines, None)
}

pub fn parse_file_to_digrams_ratio_based(
    file: &PathBuf,
    should_print_other_lines: bool,
    ratio: f32,
) -> Result<(
    NamedPath,
    NodeRegistry,
    DeterministicHashMap<usize, NodeId>,
    Vec<bool>,
)> {
    let data = bufreader_from_compressed(file)?;
    let (node_ids_by_name, number_of_paths) = parse_node_ids(data, false)?;
    let mut compressed_paths = vec![false; number_of_paths];
    let mut to_compress: DeterministicHashSet<usize> = DeterministicHashSet::default();
    let mut rng = rand::rngs::SmallRng::seed_from_u64(0);
    while to_compress.len() < (ratio * number_of_paths as f32) as usize {
        let num = rng.random_range(0..number_of_paths);
        to_compress.insert(num);
    }
    to_compress
        .into_iter()
        .for_each(|idx| compressed_paths[idx] = true);
    let data = bufreader_from_compressed(file)?;
    let (named_paths, node_ids_by_name, singleton_paths) = parse_file_content_to_digrams(
        data,
        node_ids_by_name,
        should_print_other_lines,
        Some(&compressed_paths),
    )?;
    Ok((
        named_paths,
        node_ids_by_name,
        singleton_paths,
        compressed_paths,
    ))
}

#[derive(Debug)]
pub enum DigramPath {
    Digrams(Vec<LocalizedDigram>),
    Monogram(NodeId),
}

pub fn get_digrams_from_haplotype(haplotype: &[NodeId]) -> DigramPath {
    if haplotype.len() == 1 {
        return DigramPath::Monogram(haplotype[0]);
    }

    let mut counter = 0;
    let mut digrams: Vec<LocalizedDigram> = Vec::new();
    haplotype.iter().tuple_windows().for_each(|(u, v)| {
        let address = Address::new(counter, counter + 1);
        counter += 1;
        let digram = Digram::new(*u, *v);
        let localized = LocalizedDigram::new(digram, address);
        digrams.push(localized);
    });
    DigramPath::Digrams(digrams)
}

type Path = (PathSegment, Vec<NodeId>);

fn parse_file_content_to_digrams<R: Read + Send>(
    reader: R,
    node_ids_by_name: NodeRegistry,
    should_print_other_lines: bool,
    to_compress: Option<&Vec<bool>>,
) -> Result<(NamedPath, NodeRegistry, DeterministicHashMap<usize, NodeId>)> {
    let line_reader = ByteLineReader::new(reader);
    let mut haplotypes: Vec<(PathSegment, Vec<LocalizedDigram>)> = Vec::new();
    let mut monogram_paths = DeterministicHashMap::default();
    let mut path_index = 0;
    line_reader.for_each(|line| {
        let h = path_from_line(&line, &node_ids_by_name);
        if let Some((haplotype_name, haplotype)) = h {
            if to_compress.is_none() || to_compress.expect("To-compress exists")[path_index] {
                match get_digrams_from_haplotype(&haplotype) {
                    DigramPath::Digrams(digrams) => {
                        haplotypes.push((haplotype_name, digrams));
                    }
                    DigramPath::Monogram(node) => {
                        haplotypes.push((haplotype_name, Vec::new()));
                        monogram_paths.insert(haplotypes.len() - 1, node);
                    }
                }
            }
            path_index += 1;
        } else if should_print_other_lines {
            print!("{}", str::from_utf8(&line).expect("Line is valid utf-8"));
        }
    });
    Ok((haplotypes, node_ids_by_name, monogram_paths))
}

fn parse_file_content_to_haplotypes<R: Read + Send>(
    reader: R,
    node_ids_by_name: NodeRegistry,
    should_print_other_lines: bool,
) -> Result<(Vec<Path>, NodeRegistry)> {
    let line_reader = ByteLineReader::new(reader);
    let mut haplotypes = Vec::new();
    line_reader.for_each(|line| {
        let h = path_from_line(&line, &node_ids_by_name);
        if let Some(haplotype) = h {
            haplotypes.push(haplotype);
        } else if should_print_other_lines && line[0] != b'Q' {
            print!("{}", str::from_utf8(&line).expect("Line is valid utf-8"));
        }
    });
    Ok((haplotypes, node_ids_by_name))
}

pub fn path_from_line(line: &[u8], node_ids_by_name: &NodeRegistry) -> Option<Path> {
    if line.is_empty() {
        return None;
    }
    let first_char = line[0];
    match first_char {
        b'P' => path_to_haplotype(line, node_ids_by_name),
        b'W' => walk_to_haplotype(line, node_ids_by_name),
        _ => None,
    }
}

fn path_to_haplotype(line: &[u8], node_ids_by_name: &NodeRegistry) -> Option<Path> {
    let (path_seg, buf_path_seg) = parse_path_identifier(line);
    let haplotype = parse_path_seq(buf_path_seg, node_ids_by_name);
    Some((path_seg, haplotype))
}

fn walk_to_haplotype(line: &[u8], node_ids_by_name: &NodeRegistry) -> Option<Path> {
    let (path_seg, buf_path_seg) = parse_walk_identifier(line);
    let haplotype = parse_walk_seq(buf_path_seg, node_ids_by_name);
    Some((path_seg, haplotype))
}

fn parse_node_ids<R: Read>(data: R, with_q: bool) -> Result<(NodeRegistry, usize)> {
    let mut node2id: NodeRegistry = NodeRegistry::new();

    log::info!("constructing indexes for node/edge IDs, node lengths, and P/W lines..");
    let mut buf = vec![];
    let mut data = BufReader::new(data);
    let mut number_of_paths = 0;
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'S' || (with_q && buf[0] == b'Q') {
            let mut iter = buf[2..].iter();
            let offset = iter.position(|&x| x == b'\t');
            if offset.is_none() {
                panic!("Line {} contains no tab", str::from_utf8(&buf[..]).unwrap());
            }
            let offset = offset.unwrap();
            if node2id
                .insert(buf[2..offset + 2].to_vec(), buf[0] == b'Q')
                .is_err()
            {
                println!("{}", str::from_utf8(&buf).unwrap());
                panic!(
                    "Segment with ID {} occurs multiple times in GFA",
                    str::from_utf8(&buf[2..offset + 2]).unwrap()
                )
            }
        } else if buf[0] == b'P' || buf[0] == b'W' {
            number_of_paths += 1;
        }
        buf.clear();
    }

    log::info!("found: {} nodes", node2id.len());
    Ok((node2id, number_of_paths))
}

pub type Grammar = DeterministicHashMap<UndirectedNodeId, (NodeId, NodeId)>;

fn parse_grammar<R: Read>(data: R, node_ids_by_name: &NodeRegistry) -> Result<Grammar> {
    let mut grammar = DeterministicHashMap::default();

    log::info!("constructing grammar...");
    let mut buf = vec![];
    let mut data = BufReader::new(data);

    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'Q' {
            let mut iter = buf[2..].iter();
            let offset = iter.position(|&x| x == b'\t').unwrap();
            let rule_name = &buf[2..offset + 2];
            let meta_node = node_ids_by_name.get_id(rule_name);
            let right_hand_side = parse_walk_seq(&buf[offset + 3..], node_ids_by_name);
            let right_hand_side = (right_hand_side[0], right_hand_side[1]);
            grammar.insert(meta_node, right_hand_side);
        }
        buf.clear();
    }

    log::info!("found: {} rules", grammar.len());
    Ok(grammar)
}

fn parse_path_identifier(data: &[u8]) -> (PathSegment, &[u8]) {
    let mut iter = data.iter();

    let start = iter.position(|&x| x == b'\t').unwrap() + 1;
    let offset = iter.position(|&x| x == b'\t').unwrap();
    let path_name = str::from_utf8(&data[start..start + offset]).unwrap();
    (
        PathSegment::from_str(path_name),
        &data[start + offset + 1..],
    )
}

pub fn parse_walk_identifier(data: &[u8]) -> (PathSegment, &[u8]) {
    let mut six_col: Vec<&str> = Vec::with_capacity(6);

    let mut it = data.iter();
    let mut i = 0;
    for _ in 0..6 {
        let j = it.position(|x| x == &b'\t').unwrap();
        six_col.push(str::from_utf8(&data[i..i + j]).unwrap());
        i += j + 1;
    }

    let seq_start = match six_col[4] {
        "*" => None,
        a => Some(usize::from_str(a).unwrap()),
    };

    let seq_end = match six_col[5] {
        "*" => None,
        a => Some(usize::from_str(a).unwrap()),
    };

    let path_seg = PathSegment::new(
        six_col[1].to_string(),
        six_col[2].to_string(),
        six_col[3].to_string(),
        seq_start,
        seq_end,
    );

    (path_seg, &data[i..])
}

fn parse_path_seq(data: &[u8], node_ids_by_name: &NodeRegistry) -> Vec<NodeId> {
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let haplotype = data[..end]
        .split(|&x| x == b',')
        .map(|current_node| {
            let current_node = current_node.trim_ascii();
            let orientation = Orientation::from_path(current_node[current_node.len() - 1]);
            let current_node = node_ids_by_name.get_id(&current_node[..current_node.len() - 1]);

            NodeId::new(current_node, orientation)
        })
        .collect_vec();
    log::debug!("parsing path sequences of size {} bytes..", end);
    haplotype
}

fn parse_walk_seq(data: &[u8], node_ids_by_name: &NodeRegistry) -> Vec<NodeId> {
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let haplotype = RE_WALK
        .captures_iter(&data[..end])
        .map(|m| {
            let orientation = Orientation::from_walk(m[1][0]);
            let current_node = node_ids_by_name.get_id(&m[2]);

            NodeId::new(current_node, orientation)
        })
        .collect_vec();
    log::debug!("parsing walk sequences of size {} bytes..", end);
    haplotype
}

#[cfg(test)]
pub fn parse_walk_seq_plainly(data: &[u8]) -> Vec<NodeId> {
    let haplotype = RE_WALK
        .captures_iter(data)
        .map(|m| {
            let orientation = Orientation::from_walk(m[1][0]);
            let current_node: u32 = str::from_utf8(&m[2])
                .expect("Valid utf-8")
                .parse()
                .expect("Can parse to u32");
            let current_node = UndirectedNodeId(current_node, false);

            NodeId::new(current_node, orientation)
        })
        .collect_vec();
    haplotype
}

#[cfg(test)]
pub fn get_haplotype_from_walk_string(
    text: &str,
    node_registry: &mut NodeRegistry,
) -> Vec<LocalizedDigram> {
    let data = text.as_bytes();
    let mut haplotype = Vec::new();
    let mut counter = 0;
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap_or(data.len());

    let prev_node = RE_WALK
        .captures_iter(&data[..end])
        .take(1)
        .collect::<Vec<_>>();
    let orientation = Orientation::from_walk(prev_node[0][1][0]);
    let prev_node = node_registry.get_inserted_if_not_exists(prev_node[0][2].to_vec());

    let mut prev_node = NodeId::new(prev_node, orientation);

    RE_WALK.captures_iter(&data[..end]).skip(1).for_each(|m| {
        let orientation = Orientation::from_walk(m[1][0]);
        let current_node = node_registry.get_inserted_if_not_exists(m[2].to_vec());

        let current_node = NodeId::new(current_node, orientation);

        let digram = Digram::new(prev_node, current_node);
        let address = Address::new(counter, counter + 1);
        counter += 1;
        let local_digram = LocalizedDigram::new(digram, address);
        haplotype.push(local_digram);

        prev_node = current_node;
    });
    haplotype
}

#[cfg(test)]
pub fn get_haplotypes_from_walk_strings(
    text: Vec<&str>,
    node_registry: &mut NodeRegistry,
) -> Vec<Vec<LocalizedDigram>> {
    text.into_iter()
        .map(|t| get_haplotype_from_walk_string(t, node_registry))
        .collect()
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::helpers::{utils::UndirectedNodeId, DeterministicHashMap};

    use super::*;
    use std::io::Cursor;

    fn get_node_registry() -> NodeRegistry {
        let mut dhm = DeterministicHashMap::default();
        dhm.insert("S1".bytes().collect_vec(), UndirectedNodeId::new(1));
        dhm.insert("S2".bytes().collect_vec(), UndirectedNodeId::new(2));
        dhm.insert("S3".bytes().collect_vec(), UndirectedNodeId::new(3));
        NodeRegistry::from(dhm)
    }

    #[test]
    fn test_parse_path_seq() {
        let reg = get_node_registry();
        let path_seq = "S1+,S1+,S2-,S3+\n".bytes().collect_vec();
        let haplotype = parse_path_seq(&path_seq, &reg);
        assert_eq!(haplotype.len(), 4);
    }

    #[test]
    fn test_parse_walk_seq() {
        let reg = get_node_registry();
        let walk_seq = ">S1>S1<S2>S3\n".bytes().collect_vec();
        let haplotype = parse_walk_seq(&walk_seq, &reg);
        assert_eq!(haplotype.len(), 4);
    }

    #[test]
    fn test_parse_path_identifier() {
        let data = "P\tA#0#ABC0.1:12-13\t12345".bytes().collect_vec();
        let (path_seg, remaining_data) = parse_path_identifier(&data);
        let expected = PathSegment::new(
            "A".to_string(),
            "0".to_string(),
            "ABC0.1".to_string(),
            Some(12),
            Some(13),
        );
        assert_eq!(path_seg, expected);
        assert_eq!(remaining_data.len(), 5);
    }

    #[test]
    fn test_parse_walk_identifier() {
        let data = "W\tA\t0\tABC0.1\t12\t13\t12345".bytes().collect_vec();
        let (path_seg, remaining_data) = parse_walk_identifier(&data);
        let expected = PathSegment::new(
            "A".to_string(),
            "0".to_string(),
            "ABC0.1".to_string(),
            Some(12),
            Some(13),
        );
        assert_eq!(path_seg, expected);
        assert_eq!(remaining_data.len(), 5);
    }

    #[test]
    fn test_parse_file_content_to_haplotypes() -> anyhow::Result<()> {
        let test_data =
            "P\tA#0#ABC0.1:12-13\tS1+,S1+,S2-,S3+\nW\tA\t0\tABC0.1\t12\t13\t>S1>S1<S2>S3\n";
        let cursor = Cursor::new(test_data);
        let node_ids_by_name = get_node_registry();
        let (haplotypes, _) = parse_file_content_to_haplotypes(cursor, node_ids_by_name, false)?;
        assert_eq!(haplotypes.len(), 2);
        Ok(())
    }
}
