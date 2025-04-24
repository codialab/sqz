use std::fmt;
use std::ops::{Add, AddAssign, Sub};

pub type RawNodeId = u64;

#[derive(Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct NodeId(u64);

impl NodeId {
    pub fn new(id: RawNodeId, orientation: u64) -> Self {
        Self(id << 1 ^ orientation)
    }
    pub fn flip(&self) -> NodeId {
        Self(self.0 ^ 1)
    }

    pub fn get_id(&self) -> RawNodeId {
        self.0 >> 1
    }

    pub fn get_orientation(&self) -> u64 {
        self.0 & 1
    }

    pub fn get_idx(&self) -> usize {
        self.0 as usize
    }

    pub fn from_raw(value: u64) -> Self {
        Self(value)
    }

    pub fn get_forward(&self) -> Self {
        Self(self.0 & (!1u64))
    }

    pub fn get_index_string(&self) -> String {
        format!("{}", self.0 >> 1)
    }

    pub fn is_forward(&self) -> bool {
        self.0 & 1 == 0
    }
}

impl Add<u64> for NodeId {
    type Output = Self;

    fn add(self, other: u64) -> Self {
        Self(self.0 + other)
    }
}

impl AddAssign<u64> for NodeId {
    fn add_assign(&mut self, other: u64) {
        *self = Self(self.0 + other);
    }
}

impl Sub for NodeId {
    type Output = u64;

    fn sub(self, other: Self) -> Self::Output {
        self.0 - other.0
    }
}

impl fmt::Debug for NodeId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}{}",
            if self.0 & 1 == 0 { '>' } else { '<' },
            (self.0 >> 1)
        )
    }
}

impl fmt::Display for NodeId {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}{}",
            if self.0 & 1 == 0 { '>' } else { '<' },
            (self.0 >> 1)
        )
    }
}
