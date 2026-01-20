use core::panic;

use crate::helpers::{
    DeterministicHashSet, NodeRegistry, Occurrence, RuleOccurrence, digram_occurrences::DigramOccurrences, utils::{Address, AddressNumber, CanonicalDigram, Digram, LocalizedDigram, NodeId}
};

pub type Rule = (NodeId, NodeId, NodeId, DeterministicHashSet<RuleOccurrence>);

fn get_canonized(
    u: NodeId,
    v: NodeId,
    haplotype: usize,
    a_0: AddressNumber,
    a_1: AddressNumber,
) -> (CanonicalDigram, Occurrence) {
    let digram = Digram::new(u, v);
    let address = Address::from_address_number(a_0, a_1);
    let (digram, address) = LocalizedDigram::new(digram, address).split_to_canonical();
    let occurrence = Occurrence::new(haplotype, address);
    (digram, occurrence)
}

fn build_rule(
    uv: CanonicalDigram,
    d: &mut DigramOccurrences,
    node_registry: &mut NodeRegistry,
) -> Rule {
    let q = node_registry.get_new_meta_node();
    let mut d_q = d.remove_digram(&uv);
    for occurrence in &d_q {
        let Some((left_neighbor, left_address_number)) = d.get_left_neighbor(&uv, occurrence)
        else {
            continue;
        };
        let (old_digram, old_occurrence) = get_canonized(
            left_neighbor,
            uv.get_u(),
            occurrence.get_haplotype(),
            left_address_number,
            occurrence.get_address().get_first(),
        );
        if let Err(e) = d.delete_occurrence(&old_digram, &old_occurrence) {
            panic!("ERROR: {:?}", e);
        }
        let (new_digram, new_occurrence) = get_canonized(
            left_neighbor,
            q,
            occurrence.get_haplotype(),
            left_address_number,
            occurrence.get_address().get_first(),
        );
        d.add_occurrence(&new_digram, new_occurrence);
    }
    if uv.get_u() == uv.get_v() {
        let (new_d_q, d_qq, d_qv) = Occurrence::split_self_loops(d_q);
        log::debug!("Self loop: {:?}{:?}", uv.get_u(), uv.get_v());
        log::debug!("\td_qq: {:?}", d_qq);
        log::debug!("\td_q: {:?}", new_d_q);
        log::debug!("\td_qv: {:?}", d_qv);
        d_q = new_d_q.into_iter().collect();
        d.add_digram(&Digram::new(q, q).into(), d_qq.into_iter().collect());
        let qv: CanonicalDigram = Digram::new(q, uv.get_v()).into();
        if qv.get_u() == q {
            d.add_digram(
                &qv,
                d_qv.into_iter().collect(),
            );
        } else {
            d.add_digram(
                &qv,
                d_qv.into_iter().map(|x| x.flip()).collect(),
            );
        }
    }
    for occurrence in &d_q {
        let Some((right_neighbor, right_address_number)) = d.get_right_neighbor(&uv, occurrence)
        else {
            continue;
        };
        let (old_digram, old_occurrence) = get_canonized(
            uv.get_v(),
            right_neighbor,
            occurrence.get_haplotype(),
            occurrence.get_address().get_second(),
            right_address_number,
        );
        if let Err(e) = d.delete_occurrence(&old_digram, &old_occurrence) {
            panic!("ERROR: {:?}", e);
        }
        let (new_digram, new_occurrence) = get_canonized(
            q,
            right_neighbor,
            occurrence.get_haplotype(),
            occurrence.get_address().get_first(),
            right_address_number,
        );
        d.add_occurrence(&new_digram, new_occurrence);
    }
    (q, uv.get_u(), uv.get_v(), d_q.into_iter().map(|o| o.into()).collect::<DeterministicHashSet<RuleOccurrence>>())
}

pub fn build_grammar(d: &mut DigramOccurrences, node_registry: &mut NodeRegistry) -> Vec<Rule> {
    let mut rules: Vec<Rule> = Vec::new();
    while let Some(uv) = d.get_most_frequent(2) {
        log::debug!("Building rule: {:?}", uv);
        rules.push(build_rule(uv, d, node_registry));
    }
    rules
}

#[cfg(test)]
mod tests {
    use crate::{helpers::utils::Orientation, parser::get_haplotypes_from_walk_strings};

    use super::*;

    #[test]
    fn test_loop_free_rule_creation() {
        let haplotypes = vec![">0>1>2>3", ">4>1>2>3", ">5>1>2>6"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 3);
    }

    #[test]
    fn test_loop_free_grammar_creation() {
        let haplotypes = vec![">0>1>2>3", ">4>1>2>3", ">5>1>2>6"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        println!("haplotypes: {:?}", haplotypes);
        let mut d = DigramOccurrences::from(haplotypes);
        println!("d: {:?}", d);
        let grammar = build_grammar(&mut d, &mut node_registry);
        assert_eq!(grammar.len(), 2);
    }

    #[test]
    fn test_self_loop_rule_creation() {
        let haplotypes = vec![">0>1>1>1>1>2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 2);
    }

    #[test]
    fn test_double_self_loop_rule_creation() {
        let haplotypes = vec![">0>1>1>1>1>1>1>1>1>2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 4);
        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 2);
    }

    #[test]
    fn test_double_self_loop_rule_creation_reverse() {
        let haplotypes = vec!["<0<1<1<1<1<1<1<1<1<2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Backward);
        let v = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Backward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u.flip());
        assert_eq!(rule.2, v.flip());
        assert_eq!(rule.3.len(), 4);
        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 2);
    }

    #[test]
    fn test_one_half_self_loop_rule_creation() {
        let haplotypes = vec![">0>1>1>1>1>1>1>2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 3);
        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_one_half_self_loop_rule_creation_reverse() {
        let haplotypes = vec!["<0<1<1<1<1<1<1<2"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 3);
        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_2() {
        let haplotypes = vec![">0>1>2>1>2>3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 2);
        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_2_reverse() {
        let haplotypes = vec!["<0<1<2<1<2<3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Backward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Backward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 2);
        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_3() {
        let haplotypes = vec![">0>1>2>1>3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 1);
        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, u).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u.flip());
        assert_eq!(rule.2, self_loop_meta.flip());
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_3_reverse() {
        let haplotypes = vec!["<0<1<2<1<3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Backward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Backward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 1);
        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, u).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u.flip());
        assert_eq!(rule.2, self_loop_meta.flip());
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_4() {
        let haplotypes = vec![">0>1>0>2>1>2>3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 1);
        let self_loop_meta = rule.0;
        let uv = Digram::new(v, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, v);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_pre_self_loop_case_4_reverse() {
        let haplotypes = vec!["<0<1<0<2<1<2<3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        let mut d = DigramOccurrences::from(haplotypes);
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Backward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Backward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 1);
        let self_loop_meta = rule.0;
        let uv = Digram::new(v, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, v);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);
    }

    #[test]
    fn test_loop_pre_loop_with_5_self_loops() {
        let haplotypes = vec![">0>1>2>1>2>1>2>1>2>1>2>3"];
        let mut node_registry = NodeRegistry::new();
        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        println!("haplos: {:?}", haplotypes);
        let mut d = DigramOccurrences::from(haplotypes);
        println!("=============");
        d.print_occurrences();
        println!("=============");
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule 1: {:?}", rule);
        println!("=============");
        d.print_occurrences();
        println!("=============");
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 5);

        let self_loop_meta = rule.0;
        let uv = Digram::new(self_loop_meta, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule 2: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, self_loop_meta);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 2);

        println!("=============");
        d.print_occurrences();
        println!("=============");
        
        assert_eq!(d.total_len(), 4);

        let qv = Digram::new(rule.0, self_loop_meta).into();
        let o = Occurrence::new(0, Address::new(9, 5));
        assert!(d.contains_occurrence(&qv, &o));
        
    }

    #[test]
    fn test_loop_reverse_forward_pre_loop() {
        let haplotypes = vec![">3>1<0<2>2>0<4"];
        let mut node_registry = NodeRegistry::new();

        node_registry.insert("0".as_bytes().to_vec(), false).unwrap();
        node_registry.insert("1".as_bytes().to_vec(), false).unwrap();
        node_registry.insert("2".as_bytes().to_vec(), false).unwrap();
        node_registry.insert("3".as_bytes().to_vec(), false).unwrap();
        node_registry.insert("4".as_bytes().to_vec(), false).unwrap();

        let haplotypes = get_haplotypes_from_walk_strings(haplotypes, &mut node_registry);
        println!("haplos: {:?}", haplotypes);
        let mut d = DigramOccurrences::from(haplotypes);
        println!("=============");
        d.print_occurrences();
        println!("=============");
        let u = NodeId::new(node_registry.get_id("0".as_bytes()), Orientation::Backward);
        let v = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Backward);
        let uv = Digram::new(u, v).into();
        assert!(!d.are_neighbors_equal());
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule 1: {:?}", rule);
        println!("=============");
        d.print_occurrences();
        println!("=============");
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, v);
        assert_eq!(rule.3.len(), 2);

        let x = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Backward);
        let y = NodeId::new(node_registry.get_id("2".as_bytes()), Orientation::Forward);
        let prohibited = Digram::new(x, y).into();
        assert!(!d.contains(&prohibited));

        let self_loop_meta = rule.0;
        let u = NodeId::new(node_registry.get_id("1".as_bytes()), Orientation::Forward);
        let uv = Digram::new(u, self_loop_meta).into();
        let rule = build_rule(uv, &mut d, &mut node_registry);
        println!("Rule 2: {:?}", rule);
        assert!(!d.are_neighbors_equal());
        assert_eq!(rule.1, u);
        assert_eq!(rule.2, self_loop_meta);
        assert_eq!(rule.3.len(), 1);

        println!("=============");
        d.print_occurrences();
        println!("=============");
        
        assert_eq!(d.total_len(), 3);
    }
}
