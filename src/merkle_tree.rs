use crate::field::FieldElement;
use num_bigint::BigUint;
use sha3::{Digest, Keccak256};

#[derive(Debug, Clone)]
pub struct MerkleTree {
    pub root: Vec<u8>,
    pub nodes: Vec<Vec<u8>>,
    pub leaf_hashes: Vec<Vec<u8>>,
    pub leaf_count: usize,
}

impl MerkleTree {
    pub fn new(leaves: &[FieldElement]) -> Self {
        for (i, leaf) in leaves.iter().enumerate() {
            println!("Leaf {}: {:?}", i, leaf.num);
        }

        let leaf_hashes: Vec<Vec<u8>> = leaves
            .iter()
            .map(|leaf| {
                let mut hasher = Keccak256::new();
                hasher.update(leaf.num.to_bytes_be());
                let hash = hasher.finalize().to_vec();
                println!("Leaf hash: {:?}", hash);
                hash
            })
            .collect();

        let mut nodes = leaf_hashes.clone();
        let mut current_layer = leaf_hashes.clone();
        let mut layer_num = 0;

        while current_layer.len() > 1 {
            let mut next_layer = vec![];
            for i in (0..current_layer.len()).step_by(2) {
                let left = &current_layer[i];
                let right = if i + 1 < current_layer.len() {
                    &current_layer[i + 1]
                } else {
                    left
                };
                let mut hasher = Keccak256::new();
                hasher.update(left);
                hasher.update(right);
                let parent_hash = hasher.finalize().to_vec();
                println!(
                    "Parent hash at index {}: {:?}, Left child index {}: {:?}, Right child index {}: {:?}",
                    i / 2,
                    parent_hash,
                    i,
                    left,
                    i + 1,
                    right
                );
                next_layer.push(parent_hash);
            }
            nodes.extend(next_layer.clone());
            current_layer = next_layer;
            layer_num += 1;
        }
        let root = current_layer[0].clone();
        println!("Merkle Tree Root: {:?}", root);

        MerkleTree {
            root,
            nodes,
            leaf_hashes,
            leaf_count: leaves.len(),
        }
    }

    pub fn get_proof(&self, index: usize) -> Vec<Vec<u8>> {
        let mut proof = vec![];
        let mut idx = index;
        let mut layer_start = 0;
        let mut layer_size = self.leaf_count;

        while layer_size > 1 {
            let sibling_idx = if idx % 2 == 0 { idx + 1 } else { idx - 1 };
            let sibling_idx = if sibling_idx >= layer_size {
                idx
            } else {
                sibling_idx
            };
            let sibling_hash = &self.nodes[layer_start + sibling_idx];
            println!(
                "Layer size: {}, Index: {}, Sibling index: {}, Sibling hash: {:?}",
                layer_size, idx, sibling_idx, sibling_hash
            );
            proof.push(sibling_hash.clone());
            idx /= 2;
            layer_start += layer_size;
            layer_size = (layer_size + 1) / 2;
        }
        proof
    }

    pub fn verify_proof(
        &self,
        leaf: &FieldElement,
        index: usize,
        proof: &[Vec<u8>],
    ) -> bool {
        let mut hasher = Keccak256::new();
        hasher.update(leaf.num.to_bytes_be());
        let mut computed_hash = hasher.finalize().to_vec();

        let mut idx = index;

        for sibling_hash in proof {
            let mut hasher = Keccak256::new();
            if idx % 2 == 0 {
                hasher.update(&computed_hash);
                hasher.update(sibling_hash);
            } else {
                hasher.update(sibling_hash);
                hasher.update(&computed_hash);
            }
            computed_hash = hasher.finalize().to_vec();
            idx /= 2;
        }

        println!("Computed root: {:?}", computed_hash);
        println!("Actual root: {:?}", self.root);

        computed_hash == self.root
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::FieldElement;
    use num_bigint::BigUint;

    #[test]
    fn test_new_merkle_tree() {
        let prime = BigUint::from(17u32);
        let leaves = [
            FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(7u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(10u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(15u32), prime.clone()).unwrap(),
        ];

        let merkle_tree = MerkleTree::new(&leaves);
        assert_eq!(merkle_tree.leaf_count, leaves.len());
        assert!(!merkle_tree.root.is_empty());
        assert_eq!(merkle_tree.leaf_hashes.len(), leaves.len());
    }

    #[test]
    fn test_get_merkle_proof() {
        let prime = BigUint::from(17u32);
        let leaves = [
            FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(7u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(10u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(15u32), prime.clone()).unwrap(),
        ];

        let merkle_tree = MerkleTree::new(&leaves);
        let proof = merkle_tree.get_proof(2);
        assert_eq!(proof.len(), 2);
    }

    #[test]
    fn test_verify_merkle_proof() {
        let prime = BigUint::from(17u32);
        let leaves = [
            FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(7u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(10u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(15u32), prime.clone()).unwrap(),
        ];

        let merkle_tree = MerkleTree::new(&leaves);
        let leaf = &leaves[2];
        let proof = merkle_tree.get_proof(2);
        let is_valid = merkle_tree.verify_proof(leaf, 2, &proof);
        assert!(is_valid);
    }
}
