use std::fmt::{self, Display};

use deepsize::DeepSizeOf;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct AddressNumber(pub u32);

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Address(pub AddressNumber, pub AddressNumber);

impl Ord for Address {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let self_norm = if self.0 <= self.1 {
            (self.0, self.1)
        } else {
            (self.1, self.0)
        };

        let other_norm = if other.0 <= other.1 {
            (other.0, other.1)
        } else {
            (other.1, other.0)
        };

        self_norm.cmp(&other_norm)
    }
}

impl PartialOrd for Address {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Debug for Address {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}|{}", self.0 .0, self.1 .0)
    }
}

impl Address {
    pub fn new(a: u32, b: u32) -> Self {
        Address(AddressNumber(a), AddressNumber(b))
    }

    pub fn from_address_number(a: AddressNumber, b: AddressNumber) -> Self {
        Self(a, b)
    }

    pub fn flip(&self) -> Self {
        Address(self.1, self.0)
    }

    pub fn get_first(&self) -> AddressNumber {
        self.0
    }

    pub fn get_second(&self) -> AddressNumber {
        self.1
    }

    pub fn is_forward(&self) -> bool {
        if self.0 == self.1 {
            panic!("Cannot determine direction of Address")
        }
        self.0 < self.1
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, DeepSizeOf)]
pub enum Orientation {
    Forward,
    Backward,
}

impl Orientation {
    pub fn flip(&self) -> Self {
        match self {
            Orientation::Forward => Orientation::Backward,
            Orientation::Backward => Orientation::Forward,
        }
    }

    pub fn from_path(byte: u8) -> Self {
        match byte {
            b'+' => Orientation::Forward,
            b'-' => Orientation::Backward,
            _ => panic!(
                "Orientation {} ({}) is not valid for paths",
                byte, byte as char
            ),
        }
    }

    pub fn from_walk(byte: u8) -> Self {
        match byte {
            b'>' => Orientation::Forward,
            b'<' => Orientation::Backward,
            _ => panic!(
                "Orientation {} ({}) is not valid for walks",
                byte, byte as char
            ),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, DeepSizeOf)]
pub struct UndirectedNodeId(pub u32, pub bool);

impl UndirectedNodeId {
    pub fn new(id: u32) -> Self {
        Self(id, false)
    }

    pub fn from(id: u32, is_meta_node: bool) -> Self {
        Self(id, is_meta_node)
    }

    pub fn new_meta_node(id: u32) -> Self {
        Self(id, true)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, DeepSizeOf)]
pub struct NodeId(pub UndirectedNodeId, pub Orientation);

impl NodeId {
    pub fn new(id_of_node: UndirectedNodeId, orientation: Orientation) -> Self {
        Self(id_of_node, orientation)
    }

    pub fn flip(&self) -> Self {
        Self(self.0, self.1.flip())
    }

    pub fn same_node(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    pub fn get_undirected(&self) -> UndirectedNodeId {
        self.0
    }

    pub fn get_forward(&self) -> NodeId {
        Self(self.0, Orientation::Forward)
    }

    pub fn is_forward(&self) -> bool {
        self.1 == Orientation::Forward
    }

    pub fn is_meta_node(&self) -> bool {
        self.0 .1
    }
}

impl Display for NodeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}",
            match self.1 {
                Orientation::Forward => ">",
                Orientation::Backward => "<",
            },
            self.0 .0
        )
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Digram(pub NodeId, pub NodeId);

impl Digram {
    pub fn new(u: NodeId, v: NodeId) -> Self {
        Self(u, v)
    }

    pub fn canonize(&self) -> Self {
        if self.0 > self.1 {
            Self(self.1.flip(), self.0.flip())
        } else if self.0.same_node(&self.1) {
            if self.0 .1 == Orientation::Backward {
                Self(self.1.flip(), self.0.flip())
            } else {
                self.clone()
            }
        } else {
            self.clone()
        }
    }

    pub fn flip(&self) -> Self {
        Self(self.1.flip(), self.0.flip())
    }

    pub fn is_canonical(&self) -> bool {
        if self.0 > self.1 {
            false
        } else if self.0.same_node(&self.1) {
            self.0 .1 != Orientation::Backward
        } else {
            true
        }
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct CanonicalDigram(pub NodeId, pub NodeId);

impl CanonicalDigram {
    pub fn get_u(&self) -> NodeId {
        self.0
    }

    pub fn get_v(&self) -> NodeId {
        self.1
    }

    // A digram is symmetric, i.e. both orientations are the same
    // if its nodes are the same but their direction is opposite
    pub fn is_symmetric(&self) -> bool {
        self.0 .0 == self.1 .0 && self.0 .1 != self.1 .1
    }
}

impl From<Digram> for CanonicalDigram {
    fn from(item: Digram) -> Self {
        let canonical = item.canonize();
        Self(canonical.0, canonical.1)
    }
}

impl From<CanonicalDigram> for Digram {
    fn from(val: CanonicalDigram) -> Self {
        Digram(val.0, val.1)
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct LocalizedDigram(pub Digram, pub Address);

impl LocalizedDigram {
    pub fn new(digram: Digram, address: Address) -> Self {
        Self(digram, address)
    }

    pub fn canonize(&self) -> Self {
        if self.0.is_canonical() {
            self.clone()
        } else {
            Self(self.0.canonize(), self.1.flip())
        }
    }

    #[allow(dead_code)]
    pub fn is_canonical(&self) -> bool {
        self.0.is_canonical()
    }

    pub fn split_to_canonical(self) -> (CanonicalDigram, Address) {
        let canonical = self.canonize();
        (canonical.0.into(), canonical.1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_same_node() {
        let a = NodeId::new(UndirectedNodeId::new(1), Orientation::Forward);
        let b = NodeId::new(UndirectedNodeId::new(1), Orientation::Backward);
        let c = NodeId::new(UndirectedNodeId::new(2), Orientation::Forward);

        assert!(a.same_node(&b));
        assert!(!a.same_node(&c));
    }

    #[test]
    fn test_digram_canonize() {
        let a = NodeId::new(UndirectedNodeId::new(2), Orientation::Forward);
        let b = NodeId::new(UndirectedNodeId::new(1), Orientation::Forward);
        let digram = Digram(a, b);
        let canonical = digram.canonize();
        assert_eq!(
            canonical,
            Digram(
                NodeId(UndirectedNodeId::new(1), Orientation::Backward),
                NodeId(UndirectedNodeId::new(2), Orientation::Backward)
            )
        )
    }
}
