use crate::field::FieldElement;
use crate::polynomial::Polynomial;
use crate::merkle_tree::MerkleTree;
use num_bigint::BigUint;

#[derive(Debug, Clone)]
pub struct FRIStep {
    pub evaluations: Vec<FieldElement>,
    pub merkle_tree: MerkleTree,
    pub domain: Vec<FieldElement>,
}

#[derive(Debug, Clone)]
pub struct FRICommitment {
    pub steps: Vec<FRIStep>,
    pub r_values: Vec<FieldElement>,
    pub final_poly: Polynomial,
}

pub fn fri_commit(
    initial_poly: &Polynomial,
    initial_domain: &[FieldElement],
    r_values: &[FieldElement],
) -> FRICommitment {
    let mut steps = vec![];
    let mut current_poly = initial_poly.clone();
    let mut current_domain = initial_domain.to_vec();

    for (round, r) in r_values.iter().enumerate() {
        println!("FRI Round {}", round);
        println!("Current polynomial: {:?}", current_poly);
        println!("Current domain size: {}", current_domain.len());

        // Evaluate polynomial on the current domain
        let evaluations = current_domain
            .iter()
            .map(|x| {
                let eval = current_poly.evaluate(x).unwrap();
                println!("P(x) at x = {:?} is {:?}", x.num, eval.num);
                eval
            })
            .collect::<Vec<_>>();

        // Make merkle tree with evaluation
        let merkle_tree = MerkleTree::new(&evaluations);

        steps.push(FRIStep {
            evaluations: evaluations.clone(),
            merkle_tree,
            domain: current_domain.clone(),
        });

        // Fold the polynomial that is low degeree
        current_poly = poly_folding(&current_poly, r);
        println!("Folded polynomial: {:?}", current_poly);

        // Update the doamin
        current_domain = update_domain(&current_domain);
        println!("Updated domain size: {}", current_domain.len());
    }

    FRICommitment {
        steps,
        r_values: r_values.to_vec(),
        final_poly: current_poly,
    }
}

fn poly_folding(poly: &Polynomial, r: &FieldElement) -> Polynomial {
    println!("Folding polynomial with r = {:?}", r);
    println!("Original polynomial coefficients: {:?}", poly.coefficients);

    let mut even_coeffs = vec![];
    let mut odd_coeffs = vec![];

    for (i, coeff) in poly.coefficients.iter().enumerate() {
        if i % 2 == 0 {
            even_coeffs.push(coeff.clone());
        } else {
            odd_coeffs.push(coeff.clone());
        }
    }

    let qe = Polynomial::new(even_coeffs);
    let qo = Polynomial::new(odd_coeffs);

    println!("Even polynomial qe(x): {:?}", qe.coefficients);
    println!("Odd polynomial qo(x): {:?}", qo.coefficients);

    // Folded Polynomial: P'(x) = qe(x^2) + r * qo(x^2)
    let scaled_qo = qo.poly_scale(&r).unwrap();
    println!("Scaled odd polynomial r * qo(x): {:?}", scaled_qo.coefficients);

    let new_poly = qe.poly_add(&scaled_qo).unwrap();
    println!("New polynomial after folding: {:?}", new_poly.coefficients);

    new_poly
}

fn update_domain(domain: &[FieldElement]) -> Vec<FieldElement> {
    // Square the element.
    // Domain size will be half.
    let half_size = domain.len() / 2;
    domain
        .iter()
        .take(half_size)
        .map(|x| x.pow(2u32).unwrap())
        .collect()
}

#[derive(Debug, Clone)]
pub struct FRIDecommitment {
    pub step_queries: Vec<FRIStepQuery>,
    pub final_poly: Polynomial,
}

#[derive(Debug, Clone)]
pub struct FRIStepQuery {
    pub x: FieldElement,
    pub evaluation: FieldElement,
    pub auth_path: Vec<Vec<u8>>,
    pub x_neg: FieldElement,
    pub evaluation_neg: FieldElement,
    pub auth_path_neg: Vec<Vec<u8>>,
}

pub fn fri_decommit(
    commitment: &FRICommitment,
    x: &FieldElement,
) -> FRIDecommitment {
    println!("FRI Decommitment");
    println!("Initial x: {:?}", x.num);

    let mut step_queries = vec![];
    let mut current_x = x.clone();

    for (i, step) in commitment.steps.iter().enumerate() {
        println!("Step {}", i);
        println!("Current x: {:?}", current_x.num);

        // Get evaluation of x and -x, and their merkle path.
        let index = step
            .domain
            .iter()
            .position(|d| d == &current_x)
            .unwrap();
        let evaluation = step.evaluations[index].clone();
        let auth_path = step.merkle_tree.get_proof(index);

        println!("P(x) at x = {:?} is {:?}", current_x.num, evaluation.num);

        let x_neg = current_x.negate().unwrap();
        let index_neg = step
            .domain
            .iter()
            .position(|d| d == &x_neg)
            .unwrap();
        let evaluation_neg = step.evaluations[index_neg].clone();
        let auth_path_neg = step.merkle_tree.get_proof(index_neg);

        println!("P(x) at -x = {:?} is {:?}", x_neg.num, evaluation_neg.num);

        step_queries.push(FRIStepQuery {
            x: current_x.clone(),
            evaluation,
            auth_path,
            x_neg,
            evaluation_neg,
            auth_path_neg,
        });

        // x = x^2
        current_x = current_x.pow(2u32).unwrap();
        println!("Updated x for next step: {:?}", current_x.num);
    }

    FRIDecommitment {
        step_queries,
        final_poly: commitment.final_poly.clone(),
    }
}

pub fn fri_verify(
    commitment: &FRICommitment,
    decommitment: &FRIDecommitment,
    x: &FieldElement,
    field_prime: &BigUint,
) -> bool {
    println!("FRI Verification");
    let mut current_x = x.clone();

    for (i, step_query) in decommitment.step_queries.iter().enumerate() {
        println!("Verifying step {}", i);
        println!("Current x: {:?}", current_x.num);

        let r = &commitment.r_values[i];
        println!("Random value r: {:?}", r.num);

        // Verify the merkle path.
        let step = &commitment.steps[i];

        let index = step
            .domain
            .iter()
            .position(|d| d == &current_x)
            .unwrap();

        let is_valid_proof = step
            .merkle_tree
            .verify_proof(&step_query.evaluation, index, &step_query.auth_path);

        if !is_valid_proof {
            println!("Invalid Merkle proof of x at step {}", i);
            return false;
        }

        let index_neg = step
            .domain
            .iter()
            .position(|d| d == &step_query.x_neg)
            .unwrap();

        let is_valid_proof_neg = step
            .merkle_tree
            .verify_proof(&step_query.evaluation_neg, index_neg, &step_query.auth_path_neg);

        if !is_valid_proof_neg {
            println!("Invalid Merkle proof of -x at step {}", i);
            return false;
        }

        // P_{i+1}(x^2) = [P_i(x) + P_i(-x)] / 2 + r * [P_i(x) - P_i(-x)] / (2x)
        let numerator1 = step_query.evaluation.add(&step_query.evaluation_neg).unwrap();
        let denominator1 = FieldElement::new(BigUint::from(2u32), field_prime.clone()).unwrap();
        let term1 = numerator1.div(&denominator1).unwrap();
        println!(
            "[P_i(x) + P_i(-x)] / 2 = ({:?} + {:?}) / 2 = {:?}",
            step_query.evaluation.num, step_query.evaluation_neg.num, term1.num
        );

        let numerator2 = step_query.evaluation.sub(&step_query.evaluation_neg).unwrap();
        let denominator2 = denominator1.mul(&current_x).unwrap();
        let term2 = numerator2.div(&denominator2).unwrap();
        let term2 = term2.mul(r).unwrap();
        println!(
            "r * [P_i(x) - P_i(-x)] / (2x) = {:?} * ({:?} - {:?}) / (2 * {:?}) = {:?}",
            r.num, step_query.evaluation.num, step_query.evaluation_neg.num, current_x.num, term2.num
        );

        let lhs = term1.add(&term2).unwrap();
        println!("Computed Pi+1(x^2): {:?}", lhs.num);

        // Compate with P_{i+1}
        if i + 1 < decommitment.step_queries.len() {
            let next_step_query = &decommitment.step_queries[i + 1];
            let expected = next_step_query.evaluation.clone();
            println!("Expected Pi+1(x^2): {:?}", expected.num);

            if lhs != expected {
                println!("Verification failed at step {}", i);
                return false;
            } else {
                println!("Step {} verification succeeded", i);
            }
        } else {
            // At final step, compare with the final_poly.
            let x_squared = current_x.pow(2u32).unwrap();
            let final_eval = decommitment.final_poly.evaluate(&x_squared).unwrap();
            println!("Final polynomial evaluation {:?}", final_eval.num);

            if lhs != final_eval {
                println!("Final verification failed");
                return false;
            } else {
                println!("Final verification succeeded");
            }
        }

        // x^2
        current_x = current_x.pow(2u32).unwrap();
        println!("Updated x for next step: {:?}", current_x.num);
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::FieldElement;
    use crate::polynomial::Polynomial;
    use num_bigint::BigUint;

    #[test]
    fn test_fri() {
        let prime = BigUint::from(17u32);

        // P0(x) = 15x^3 + 10x^2 + 15x + 1
        let coefficients = vec![
            FieldElement::new(1u32.into(), prime.clone()).unwrap(),
            FieldElement::new(15u32.into(), prime.clone()).unwrap(),
            // FieldElement::new(10u32.into(), prime.clone()).unwrap(),
            FieldElement::new(0u32.into(), prime.clone()).unwrap(),
            FieldElement::new(15u32.into(), prime.clone()).unwrap(),
        ];
        let p0 = Polynomial::new(coefficients);
        println!("Initial polynomial P0(x): {:?}", p0.coefficients);

        // generat for low degree extention, g = 9.
        let g = FieldElement::new(9u32.into(), prime.clone()).unwrap();

        // D0 = 3 * g^i mod 17 (i = 0..7).
        let mut d0 = vec![];
        let three = FieldElement::new(3u32.into(), prime.clone()).unwrap();
        for i in 0..8 {
            let g_pow = g.pow(i as u32).unwrap();
            let x = three.mul(&g_pow).unwrap();
            d0.push(x);
        }
        println!("Initial domain D0: {:?}", d0.iter().map(|x| x.num.clone()).collect::<Vec<_>>());

        // r_values = [4, 3]
        let r_values = vec![
            FieldElement::new(BigUint::from(4u32), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap(),
        ];
        println!("Random values r: {:?}", r_values.iter().map(|r| r.num.clone()).collect::<Vec<_>>());

        let commitment = fri_commit(&p0, &d0, &r_values);

        // P0(x), x = 10
        let x = FieldElement::new(11u32.into(), prime.clone()).unwrap();
        println!("Selected x: {:?}", x.num);
        let decommitment = fri_decommit(&commitment, &x);

        let is_valid = fri_verify(&commitment, &decommitment, &x, &prime);
        assert!(is_valid, "FRI proof is invalid");
        println!("FRI proof is valid");
    }
}
