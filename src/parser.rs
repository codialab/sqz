use anyhow::{anyhow, Result};
use flate2::read::MultiGzDecoder;
use lazy_static::lazy_static;
use regex::bytes::Regex;
use std::{collections::HashMap, io::BufReader, str::FromStr};

use crate::helpers::{
    utils::Address, utils::Digram, utils::LocalizedDigram, utils::NodeId, utils::Orientation,
    utils::UndirectedNodeId, Haplotype, NodeRegistry, PathSegment,
};
use std::{
    collections::HashSet,
    io::{self, BufRead, Read},
    path::PathBuf,
};

lazy_static! {
    static ref RE_WALK: Regex = Regex::new(r"([><])([!-;=?-~]+)").unwrap();
}

struct ByteLineReader<R: io::Read> {
    data: io::BufReader<R>,
}

impl<R: io::Read> ByteLineReader<R> {
    fn new(data: R) -> Self {
        Self {
            data: BufReader::new(data),
        }
    }
}

impl<R: io::Read> Iterator for ByteLineReader<R> {
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

fn bufreader_from_compressed(file: &PathBuf) -> Result<io::BufReader<Box<dyn Read>>> {
    let f = std::fs::File::open(file)?;
    if let Some(name) = file.file_name() {
        let reader: Box<dyn Read> = if name.to_str().unwrap().ends_with(".gz") {
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

pub fn parse_file_to_haplotypes(
    file: &PathBuf,
    should_print_other_lines: bool,
) -> Result<(Vec<Haplotype>, NodeRegistry, HashMap<usize, NodeId>)> {
    let data = bufreader_from_compressed(file)?;
    let node_ids_by_name = parse_node_ids(data, false)?;
    let data = bufreader_from_compressed(file)?;
    parse_file_content_to_haplotypes(data, node_ids_by_name, should_print_other_lines)
}

fn parse_file_content_to_haplotypes<R: Read>(
    reader: R,
    node_ids_by_name: NodeRegistry,
    should_print_other_lines: bool,
) -> Result<(Vec<Haplotype>, NodeRegistry, HashMap<usize, NodeId>)> {
    let line_reader = ByteLineReader::new(reader);
    let mut counter = 0;
    let mut haplotypes = Vec::new();
    let mut single_node_haplotypes = HashMap::new();
    line_reader.for_each(|line| {
        let h = line_to_haplotype(&line, &node_ids_by_name);
        if let Some((haplotype, first_node)) = h {
            if haplotype.1.is_empty() {
                single_node_haplotypes.insert(counter, first_node);
            }
            haplotypes.push(haplotype);
            counter += 1;
        } else if should_print_other_lines {
            print!("{}", str::from_utf8(&line).expect("Line is valid utf-8"));
        }
    });
    Ok((haplotypes, node_ids_by_name, single_node_haplotypes))
}

fn line_to_haplotype(line: &[u8], node_ids_by_name: &NodeRegistry) -> Option<(Haplotype, NodeId)> {
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

fn path_to_haplotype(line: &[u8], node_ids_by_name: &NodeRegistry) -> Option<(Haplotype, NodeId)> {
    let (path_seg, buf_path_seg) = parse_path_identifier(line);
    let (haplotype, first_node) = parse_path_seq(buf_path_seg, node_ids_by_name);
    Some(((path_seg, haplotype), first_node))
}

fn walk_to_haplotype(line: &[u8], node_ids_by_name: &NodeRegistry) -> Option<(Haplotype, NodeId)> {
    let (path_seg, buf_path_seg) = parse_walk_identifier(line);
    let (haplotype, first_node) = parse_walk_seq(buf_path_seg, node_ids_by_name);
    Some(((path_seg, haplotype), first_node))
}

fn parse_node_ids<R: Read>(data: R, with_q: bool) -> Result<NodeRegistry> {
    let mut node2id: NodeRegistry = NodeRegistry::new();

    log::info!("constructing indexes for node/edge IDs, node lengths, and P/W lines..");
    let mut buf = vec![];
    let mut data = BufReader::new(data);
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'S' || (with_q && buf[0] == b'Q') {
            let mut iter = buf[2..].iter();
            let offset = iter.position(|&x| x == b'\t').unwrap();
            if node2id.insert(buf[2..offset + 2].to_vec()).is_err() {
                println!("{}", str::from_utf8(&buf).unwrap());
                panic!(
                    "Segment with ID {} occurs multiple times in GFA",
                    str::from_utf8(&buf[2..offset + 2]).unwrap()
                )
            }
        }
        buf.clear();
    }

    log::info!("found: {} nodes", node2id.len());
    Ok(node2id)
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

fn parse_path_seq(data: &[u8], node_ids_by_name: &NodeRegistry) -> (Vec<LocalizedDigram>, NodeId) {
    let mut haplotype = Vec::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<UndirectedNodeId> = HashSet::new();
    let mut curr_counter = 0;
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let prev_node = data[..end]
        .split(|&x| x == b',')
        .take(1)
        .collect::<Vec<_>>()[0];
    let prev_node = prev_node.trim_ascii();
    let orientation = Orientation::from_path(prev_node[prev_node.len() - 1]);
    let prev_node = node_ids_by_name.get_id(&prev_node[..prev_node.len() - 1]);

    // Set the counter before setting the orientation to not take the orientation into account
    let mut prev_node = NodeId::new(prev_node, orientation);

    let first_node = prev_node;

    nodes_visited_curr.insert(prev_node.get_undirected());
    data[..end]
        .split(|&x| x == b',')
        .skip(1)
        .for_each(|current_node| {
            let current_node = current_node.trim_ascii();
            let orientation = Orientation::from_path(current_node[current_node.len() - 1]);
            let current_node = node_ids_by_name.get_id(&current_node[..current_node.len() - 1]);

            let current_node = NodeId::new(current_node, orientation);

            prev_counter = curr_counter;
            if nodes_visited_curr.contains(&current_node.get_undirected()) {
                nodes_visited_curr.clear();
            }
            curr_counter += 1;
            nodes_visited_curr.insert(current_node.get_undirected());

            let digram = Digram::new(prev_node, current_node);
            let address = Address::new(prev_counter, curr_counter);
            let local_digram = LocalizedDigram::new(digram, address);
            haplotype.push(local_digram);

            prev_node = current_node;
        });
    log::debug!("parsing path sequences of size {} bytes..", end);
    (haplotype, first_node)
}

fn parse_walk_seq(data: &[u8], node_ids_by_name: &NodeRegistry) -> (Vec<LocalizedDigram>, NodeId) {
    let mut haplotype = Vec::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<UndirectedNodeId> = HashSet::new();
    let mut curr_counter = 0;
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    let prev_node = RE_WALK
        .captures_iter(&data[..end])
        .take(1)
        .collect::<Vec<_>>();
    let orientation = Orientation::from_walk(prev_node[0][1][0]);
    let prev_node = node_ids_by_name.get_id(&prev_node[0][2]);

    let mut prev_node = NodeId::new(prev_node, orientation);

    let first_node = prev_node;

    nodes_visited_curr.insert(prev_node.get_undirected());

    RE_WALK.captures_iter(&data[..end]).skip(1).for_each(|m| {
        let orientation = Orientation::from_walk(m[1][0]);
        let current_node = node_ids_by_name.get_id(&m[2]);

        let current_node = NodeId::new(current_node, orientation);

        prev_counter = curr_counter;
        if nodes_visited_curr.contains(&current_node.get_undirected()) {
            curr_counter += 1;
            nodes_visited_curr.clear();
        }
        nodes_visited_curr.insert(current_node.get_undirected());

        let digram = Digram::new(prev_node, current_node);
        let address = Address::new(prev_counter, curr_counter);
        let local_digram = LocalizedDigram::new(digram, address);
        haplotype.push(local_digram);

        prev_node = current_node;
    });
    log::debug!("parsing walk sequences of size {} bytes..", end);
    (haplotype, first_node)
}

#[cfg(test)]
pub fn get_haplotype_from_walk_string(
    text: &str,
    node_registry: &mut NodeRegistry,
) -> Vec<LocalizedDigram> {
    let data = text.as_bytes();
    let mut haplotype = Vec::new();
    let mut prev_counter = 0;
    let mut nodes_visited_curr: HashSet<UndirectedNodeId> = HashSet::new();
    let mut curr_counter = 0;
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

    nodes_visited_curr.insert(prev_node.get_undirected());

    RE_WALK.captures_iter(&data[..end]).skip(1).for_each(|m| {
        let orientation = Orientation::from_walk(m[1][0]);
        let current_node = node_registry.get_inserted_if_not_exists(m[2].to_vec());

        let current_node = NodeId::new(current_node, orientation);

        prev_counter = curr_counter;
        if nodes_visited_curr.contains(&current_node.get_undirected()) {
            curr_counter += 1;
            nodes_visited_curr.clear();
        }
        nodes_visited_curr.insert(current_node.get_undirected());

        let digram = Digram::new(prev_node, current_node);
        let address = Address::new(prev_counter, curr_counter);
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

    use crate::helpers::DeterministicHashMap;

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
        let (haplotype, _) = parse_path_seq(&path_seq, &reg);
        assert_eq!(haplotype.len(), 3);
    }

    #[test]
    fn test_parse_walk_seq() {
        let reg = get_node_registry();
        let walk_seq = ">S1>S1<S2>S3\n".bytes().collect_vec();
        let (haplotype, _) = parse_walk_seq(&walk_seq, &reg);
        assert_eq!(haplotype.len(), 3);
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
        let (haplotypes, _, _) = parse_file_content_to_haplotypes(cursor, node_ids_by_name, false)?;
        assert_eq!(haplotypes.len(), 2);
        Ok(())
    }
}
