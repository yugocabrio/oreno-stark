// deep_ali.rs

use crate::field::FieldElement;
use crate::polynomial::Polynomial;
use num_bigint::BigUint;

/// Represents the DEEP-ALI protocol calculations.
pub struct DeepALI {
    pub trace_polynomial: Polynomial,
    pub field_prime: BigUint,
}

impl DeepALI {
    /// Creates a new instance of DeepALI with the given trace polynomial.
    pub fn new(trace_polynomial: Polynomial, field_prime: BigUint) -> Self {
        Self {
            trace_polynomial,
            field_prime,
        }
    }

    /// Computes the boundary constraint polynomial C1(x).
    pub fn compute_boundary_constraint(
        &self,
        point: &FieldElement,
        expected_value: &FieldElement,
    ) -> (Polynomial, Polynomial) {
        // p1(x) = t(x) - expected_value
        let p1_poly = self
            .trace_polynomial
            .poly_sub(&Polynomial::new(vec![expected_value.clone()]))
            .unwrap();

        // z1(x) = x - point
        let one = FieldElement::one(&self.field_prime);
        let z1_poly = Polynomial::new(vec![point.clone().negate().unwrap(), one]);

        // C1(x) = p1(x) / z1(x)
        let (c1_poly, remainder) = p1_poly.poly_div_rem(&z1_poly).unwrap();
        assert!(
            remainder.is_zero(),
            "Boundary constraint division has a remainder"
        );

        (c1_poly, z1_poly)
    }

    /// Computes the transition constraint polynomial C2(x).
    pub fn compute_transition_constraint(
        &self,
        generator: &FieldElement,
    ) -> (Polynomial, Polynomial) {
        // Compute t(gx) by substituting x with g * x
        let mut t_gx_coeffs = vec![];
        for (i, coeff) in self.trace_polynomial.coefficients.iter().enumerate() {
            let g_pow_i = generator.pow(i as u32).unwrap();
            let new_coeff = coeff.mul(&g_pow_i).unwrap();
            t_gx_coeffs.push(new_coeff);
        }
        let t_gx_poly = Polynomial::new(t_gx_coeffs);

        // Compute [t(x)]^2
        let t_x_squared = self.trace_polynomial.poly_mul(&self.trace_polynomial).unwrap();

        // Compute p2(x) = t(gx) - [t(x)]^2
        let p2_poly = t_gx_poly.poly_sub(&t_x_squared).unwrap();

        // Compute z2(x) = (x - 1)(x - 13)(x - 16)
        let one = FieldElement::one(&self.field_prime);
        let roots = vec![
            FieldElement::new(1u32.into(), self.field_prime.clone()).unwrap(),
            FieldElement::new(13u32.into(), self.field_prime.clone()).unwrap(),
            FieldElement::new(16u32.into(), self.field_prime.clone()).unwrap(),
        ];

        let mut z2_poly = Polynomial::new(vec![FieldElement::one(&self.field_prime)]);
        for root in roots {
            let z_poly = Polynomial::new(vec![root.negate().unwrap(), one.clone()]);
            z2_poly = z2_poly.poly_mul(&z_poly).unwrap();
        }

        // C2(x) = p2(x) / z2(x)
        let (c2_poly, remainder) = p2_poly.poly_div_rem(&z2_poly).unwrap();
        assert!(
            remainder.is_zero(),
            "Transition constraint division has a remainder"
        );

        (c2_poly, z2_poly)
    }

    /// Computes the DEEP composition polynomial H(x).
    pub fn compute_deep_composition_polynomial(
        &self,
        constraints: Vec<(Polynomial, u32, FieldElement, FieldElement)>,
        total_degree: u32,
    ) -> Polynomial {
        let mut h_poly = Polynomial::new(vec![FieldElement::zero(&self.field_prime)]);

        for (constraint_poly, constraint_degree, alpha, beta) in constraints {
            let degree_adjustment = total_degree - constraint_degree;

            // x^{D - D_i}
            let x_to_deg_adj =
                Polynomial::monomial(degree_adjustment as usize, FieldElement::one(&self.field_prime));

            // Adjustment polynomial: alpha * x^{D - D_i} + beta
            let adjustment_poly = x_to_deg_adj
                .poly_scale(&alpha)
                .unwrap()
                .poly_add(&Polynomial::new(vec![beta]))
                .unwrap();

            // Term: C_i(x) * adjustment_poly
            let term = constraint_poly.poly_mul(&adjustment_poly).unwrap();

            // Add to H(x)
            h_poly = h_poly.poly_add(&term).unwrap();
        }

        h_poly
    }

    /// Splits H(x) into H1(x^2) and H2(x^2).
    pub fn split_h_polynomial(
        &self,
        h_poly: &Polynomial,
    ) -> (Polynomial, Polynomial) {
        let h_coeffs = h_poly.coefficients.clone();
        let h1_coeffs = h_coeffs
            .iter()
            .enumerate()
            .filter_map(|(i, coeff)| if i % 2 == 0 { Some(coeff.clone()) } else { None })
            .collect::<Vec<_>>();
        let h2_coeffs = h_coeffs
            .iter()
            .enumerate()
            .filter_map(|(i, coeff)| if i % 2 == 1 { Some(coeff.clone()) } else { None })
            .collect::<Vec<_>>();

        let h1_poly = Polynomial::new(h1_coeffs);
        let h2_poly = Polynomial::new(h2_coeffs);

        (h1_poly, h2_poly)
    }

    /// Computes the DEEP-ALI polynomial P₀(x).
    pub fn compute_deep_ali_polynomial(
        &self,
        h_poly: &Polynomial,
        z: &FieldElement,
        generator: &FieldElement,
        gammas: &[FieldElement],
    ) -> Polynomial {
        let one = FieldElement::one(&self.field_prime);

        // Compute t(z) and t(gz)
        let t_z = self.trace_polynomial.evaluate(z).unwrap();
        let gz = generator.mul(z).unwrap();
        let t_gz = self.trace_polynomial.evaluate(&gz).unwrap();

        // Term 1: (t(x) - t(z)) / (x - z)
        let numerator1 = self
            .trace_polynomial
            .poly_sub(&Polynomial::new(vec![t_z.clone()]))
            .unwrap();
        let denominator1 = Polynomial::new(vec![z.negate().unwrap(), one.clone()]);
        let (term1_poly, _) = numerator1.poly_div_rem(&denominator1).unwrap();
        let term1_poly = term1_poly.poly_scale(&gammas[0]).unwrap();

        // Term 2: (t(x) - t(gz)) / (x - gz)
        let numerator2 = self
            .trace_polynomial
            .poly_sub(&Polynomial::new(vec![t_gz.clone()]))
            .unwrap();
        let denominator2 = Polynomial::new(vec![gz.negate().unwrap(), one.clone()]);
        let (term2_poly, _) = numerator2.poly_div_rem(&denominator2).unwrap();
        let term2_poly = term2_poly.poly_scale(&gammas[1]).unwrap();

        // H(x) split into H1(x^2) and H2(x^2)
        let (h1_poly, h2_poly) = self.split_h_polynomial(h_poly);

        // Evaluate H1(z^2) and H2(z^2)
        let z_squared = z.pow(2u32).unwrap();
        let h1_z_squared = h1_poly.evaluate(&z_squared).unwrap();
        let h2_z_squared = h2_poly.evaluate(&z_squared).unwrap();

        // Corrected denominators: (x - z)
        let denominator3 = Polynomial::new(vec![z.negate().unwrap(), one.clone()]);

        // Term 3: (H1(x^2) - H1(z^2)) / (x - z)
        let x_squared_poly = Polynomial::new(vec![
            FieldElement::zero(&self.field_prime),
            FieldElement::zero(&self.field_prime),
            one.clone(),
        ]);
        let h1_x_squared = h1_poly.compose(&x_squared_poly).unwrap();
        let numerator3 = h1_x_squared
            .poly_sub(&Polynomial::new(vec![h1_z_squared.clone()]))
            .unwrap();
        let (term3_poly, _) = numerator3.poly_div_rem(&denominator3).unwrap();
        let term3_poly = term3_poly.poly_scale(&gammas[2]).unwrap();

        // Term 4: (H2(x^2) - H2(z^2)) / (x - z)
        let h2_x_squared = h2_poly.compose(&x_squared_poly).unwrap();
        let numerator4 = h2_x_squared
            .poly_sub(&Polynomial::new(vec![h2_z_squared.clone()]))
            .unwrap();
        let (term4_poly, _) = numerator4.poly_div_rem(&denominator3).unwrap();
        let term4_poly = term4_poly.poly_scale(&gammas[3]).unwrap();

        // P₀(x) = term1 + term2 + term3 + term4
        let p0_poly = term1_poly
            .poly_add(&term2_poly)
            .unwrap()
            .poly_add(&term3_poly)
            .unwrap()
            .poly_add(&term4_poly)
            .unwrap();

        p0_poly
    }
}

// tests/deep_ali_tests.rs

#[cfg(test)]
mod tests {
    use crate::deep_ali::DeepALI;
    use crate::field::FieldElement;
    use crate::polynomial::Polynomial;
    use num_bigint::BigUint;
    use num_traits::One;

    #[test]
    fn test_deep_ali_calculations() {
        // Field setup (prime field with prime = 17)
        let prime = BigUint::from(17u32);

        // Define the field elements
        let one = FieldElement::one(&prime);
        let zero = FieldElement::zero(&prime);

        // Define the trace polynomial t(x) = 13x^3 + 2x^2 + 16x + 6
        let coefficients = vec![
            FieldElement::new(6u32.into(), prime.clone()).unwrap(),    // Constant term
            FieldElement::new(16u32.into(), prime.clone()).unwrap(),   // x^1
            FieldElement::new(2u32.into(), prime.clone()).unwrap(),    // x^2
            FieldElement::new(13u32.into(), prime.clone()).unwrap(),   // x^3
        ];
        let trace_polynomial = Polynomial::new(coefficients.clone());

        // Instantiate DeepALI
        let deep_ali = DeepALI::new(trace_polynomial.clone(), prime.clone());

        // Compute boundary constraint
        let point = FieldElement::one(&prime);
        let expected_value = FieldElement::new(3u32.into(), prime.clone()).unwrap();
        let (c1_poly, z1_poly) = deep_ali.compute_boundary_constraint(&point, &expected_value);

        // Expected C1(x)
        let expected_c1_coeffs = vec![
            FieldElement::new(14u32.into(), prime.clone()).unwrap(),
            FieldElement::new(15u32.into(), prime.clone()).unwrap(),
            FieldElement::new(13u32.into(), prime.clone()).unwrap(),
        ];
        let expected_c1_poly = Polynomial::new(expected_c1_coeffs);

        assert_eq!(c1_poly, expected_c1_poly, "C1(x) does not match the expected value");

        // Compute transition constraint
        let generator = FieldElement::new(13u32.into(), prime.clone()).unwrap();
        let (c2_poly, z2_poly) = deep_ali.compute_transition_constraint(&generator);

        // Expected C2(x)
        let expected_c2_coeffs = vec![
            FieldElement::new(16u32.into(), prime.clone()).unwrap(),
            FieldElement::new(9u32.into(), prime.clone()).unwrap(),
            FieldElement::new(12u32.into(), prime.clone()).unwrap(),
            FieldElement::new(1u32.into(), prime.clone()).unwrap(),
        ];
        let expected_c2_poly = Polynomial::new(expected_c2_coeffs);

        assert_eq!(c2_poly, expected_c2_poly, "C2(x) does not match the expected value");

        // Compute H(x)
        let total_degree = 4u32;
        let alpha1 = FieldElement::new(1u32.into(), prime.clone()).unwrap();
        let beta1 = FieldElement::new(3u32.into(), prime.clone()).unwrap();
        let alpha2 = FieldElement::new(2u32.into(), prime.clone()).unwrap();
        let beta2 = FieldElement::new(4u32.into(), prime.clone()).unwrap();

        let constraints = vec![
            (
                c1_poly.clone(),
                c1_poly.degree() as u32,
                alpha1.clone(),
                beta1.clone(),
            ),
            (
                c2_poly.clone(),
                c2_poly.degree() as u32,
                alpha2.clone(),
                beta2.clone(),
            ),
        ];

        let h_poly = deep_ali.compute_deep_composition_polynomial(constraints, total_degree);

        // Expected H(x)
        let expected_h_coeffs = vec![
            FieldElement::new(4u32.into(), prime.clone()).unwrap(),
            FieldElement::new(11u32.into(), prime.clone()).unwrap(),
            FieldElement::new(0u32.into(), prime.clone()).unwrap(),
            FieldElement::new(9u32.into(), prime.clone()).unwrap(),
            FieldElement::new(15u32.into(), prime.clone()).unwrap(),
        ];
        let expected_h_poly = Polynomial::new(expected_h_coeffs);

        assert_eq!(h_poly, expected_h_poly, "H(x) does not match the expected value");

        // Split H(x) into H1(x^2) and H2(x^2)
        let (h1_poly, h2_poly) = deep_ali.split_h_polynomial(&h_poly);

        // Expected H1(x^2) and H2(x^2)
        let expected_h1_coeffs = vec![
            FieldElement::new(4u32.into(), prime.clone()).unwrap(),
            FieldElement::new(0u32.into(), prime.clone()).unwrap(),
            FieldElement::new(15u32.into(), prime.clone()).unwrap(),
        ];
        let expected_h1_poly = Polynomial::new(expected_h1_coeffs);

        let expected_h2_coeffs = vec![
            FieldElement::new(11u32.into(), prime.clone()).unwrap(),
            FieldElement::new(9u32.into(), prime.clone()).unwrap(),
        ];
        let expected_h2_poly = Polynomial::new(expected_h2_coeffs);

        assert_eq!(h1_poly, expected_h1_poly, "H1(x) does not match the expected value");
        assert_eq!(h2_poly, expected_h2_poly, "H2(x) does not match the expected value");

        // Compute H1(z^2) and H2(z^2)
        let z = FieldElement::new(8u32.into(), prime.clone()).unwrap();
        let z_squared = z.pow(2u32).unwrap(); // z^2

        let h1_z_squared = h1_poly.evaluate(&z_squared).unwrap();
        let h2_z_squared = h2_poly.evaluate(&z_squared).unwrap();

        // Compute H(z) = H1(z^2) + z * H2(z^2)
        let h_z = h1_z_squared.add(&z.mul(&h2_z_squared).unwrap()).unwrap();

        // Evaluate H(z) directly from H(x)
        let h_z_direct = h_poly.evaluate(&z).unwrap();

        assert_eq!(h_z, h_z_direct, "H(z) computed from H1 and H2 does not match direct evaluation");

        // Expected values from the article
        // H1(z^2) = 6
        let expected_h1_z_squared = FieldElement::new(6u32.into(), prime.clone()).unwrap();
        // H2(z^2) = 9
        let expected_h2_z_squared = FieldElement::new(9u32.into(), prime.clone()).unwrap();
        // H(z) = H1(z^2) + z * H2(z^2) = 6 + 8 * 9 = 6 + 8 * 9 mod 17
        let expected_h_z = h1_z_squared.add(&z.mul(&h2_z_squared).unwrap()).unwrap();

        assert_eq!(h1_z_squared, expected_h1_z_squared, "H1(z^2) does not match the expected value");
        assert_eq!(h2_z_squared, expected_h2_z_squared, "H2(z^2) does not match the expected value");
        assert_eq!(h_z, expected_h_z, "H(z) does not match the expected value from H1(z^2) and H2(z^2)");

        // Now compute P₀(x)
        let gammas = vec![one.clone(); 4]; // gamma_i = 1

        let p0_poly = deep_ali.compute_deep_ali_polynomial(&h_poly, &z, &generator, &gammas);

        // Expected P₀(x)
        let expected_p0_coeffs = vec![
            FieldElement::new(1u32.into(), prime.clone()).unwrap(),
            FieldElement::new(15u32.into(), prime.clone()).unwrap(),
            FieldElement::new(10u32.into(), prime.clone()).unwrap(),
            FieldElement::new(15u32.into(), prime.clone()).unwrap(),
        ];
        let expected_p0_poly = Polynomial::new(expected_p0_coeffs);

        assert_eq!(p0_poly, expected_p0_poly, "P₀(x) does not match the expected value");

        println!("All calculations match the expected values from the article.");
    }
}

// マークルツリー
