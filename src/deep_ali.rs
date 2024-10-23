// TODO
// Write fnctions that converts Trace Polynomial → DEEP-ALI.

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::{FromPrimitive, Zero, One};
    use crate::field::FieldElement;
    use crate::polynomial::Polynomial;

    // DEEP ALI calculations by numerically.
    // The values in the test functions follow the examples in the following article.
    // https://blog.lambdaclass.com/diving-deep-fri/
    #[test]
    fn test_deep_ali() {
        let prime = BigUint::from_u32(17u32).unwrap();

        let generator_value = BigUint::from_u32(3u32).unwrap();
        let generator_element = FieldElement::new(generator_value.clone(), prime.clone()).unwrap();
        let generator = generator_element.set_generator(prime.clone()).unwrap();

        let n = 4u32;
        let n_biguint = BigUint::from_u32(n).unwrap();
        let omega = generator.nth_generator(n_biguint.clone()).unwrap();

        println!("4th omega: {:?}", omega);

        // t(x) = 13x^3 + 2x^2 + 16x + 6
        let t_coefficients = vec![
            FieldElement::new(BigUint::from_u32(6u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(2u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap(),
        ];
        let t_polynomial = Polynomial::new(t_coefficients.clone());

        // Boundary constraints
        // p1(x) = t(x) - 3
        let three = FieldElement::new(BigUint::from_u32(3u32).unwrap(), prime.clone()).unwrap();
        let mut t_minus_three_coeffs = t_coefficients.clone();
        t_minus_three_coeffs[0] = t_minus_three_coeffs[0].sub(&three).unwrap();
        let p1_polynomial = Polynomial::new(t_minus_three_coeffs);

        // z1(x) = x - 1
        let one = FieldElement::new(BigUint::from_u32(1u32).unwrap(), prime.clone()).unwrap();
        let z1_polynomial = Polynomial::new(vec![
            one.clone().negate().unwrap(),
            FieldElement::new(BigUint::one(), prime.clone()).unwrap(),
        ]);

        // c1(x) = p1(x)/z1(x)
        let (c1_polynomial, remainder) = p1_polynomial.poly_div_rem(&z1_polynomial).unwrap();
        assert!(remainder.coefficients.is_empty() || remainder.coefficients.iter().all(|c| c.num.is_zero()));

        println!("C1(x): {:?}", c1_polynomial);

        // Transition constraint 
        // p2(x) = t(gx) - t(x)^2

        // t(gx)
        let gx_coeffs = t_polynomial.coefficients.iter().enumerate().map(|(i, coeff)| {
            let omega_pow_i = omega.pow(i as u32).unwrap();
            coeff.mul(&omega_pow_i).unwrap()
        }).collect::<Vec<_>>();
        let t_gx_polynomial = Polynomial::new(gx_coeffs);
        
        println!("t_gx_polynomial : {:?}", t_gx_polynomial);

        // t(x)^2
        let t_x_square = t_polynomial.poly_mul(&t_polynomial).unwrap();

        println!("t_x_square : {:?}", t_x_square);

        // p2(x) = t(gx) - t(x)^2
        let p2_polynomial = t_gx_polynomial.poly_sub(&t_x_square).unwrap();
        println!("p2_polynomial : {:?}", p2_polynomial);

    
        // z2(x) = (x - 1)(x - 13)(x - 16)
        let thirteen = FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap();
        let sixteen = FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap();

        let x_minus_one = Polynomial::new(vec![
            one.clone().negate().unwrap(),
            FieldElement::new(BigUint::one(), prime.clone()).unwrap(),
        ]);

        let x_minus_thirteen = Polynomial::new(vec![
            thirteen.clone().negate().unwrap(),
            FieldElement::new(BigUint::one(), prime.clone()).unwrap(),
        ]);

        let x_minus_sixteen = Polynomial::new(vec![
            sixteen.clone().negate().unwrap(),
            FieldElement::new(BigUint::one(), prime.clone()).unwrap(),
        ]);

        let z2_polynomial = x_minus_one.poly_mul(&x_minus_thirteen).unwrap().poly_mul(&x_minus_sixteen).unwrap();

        println!("z2_polynomial : {:?}", z2_polynomial);

        // c2(x) = p2(x) / x2(x)
        let (c2_polynomial, remainder_p2) = p2_polynomial.poly_div_rem(&z2_polynomial).unwrap();
        assert!(remainder_p2.coefficients.is_empty() || remainder_p2.coefficients.iter().all(|c| c.num.is_zero()));

        println!("C2(x): {:?}", c2_polynomial);

        // The (constraint) composition polynomial
        // α1, β1, α2, β2
        let alpha1 = FieldElement::new(BigUint::from_u32(1u32).unwrap(), prime.clone()).unwrap();
        let beta1 = FieldElement::new(BigUint::from_u32(3u32).unwrap(), prime.clone()).unwrap();
        let alpha2 = FieldElement::new(BigUint::from_u32(2u32).unwrap(), prime.clone()).unwrap();
        let beta2 = FieldElement::new(BigUint::from_u32(4u32).unwrap(), prime.clone()).unwrap();

        let d = 4u32;
        let d1 = 2u32;

        // C1(x) * (α1*x^(D-D1) + β1)
        // x^(D-D1)
        let exponent = d - d1;
        let mut x_to_d_coeffs = vec![FieldElement::new(BigUint::zero(), prime.clone()).unwrap(); exponent as usize];
        x_to_d_coeffs.push(FieldElement::new(BigUint::one(), prime.clone()).unwrap());
        let x_to_d = Polynomial::new(x_to_d_coeffs);

        let alpha1_xd = x_to_d.poly_scale(&alpha1).unwrap();
        let beta1_poly = Polynomial::new(vec![beta1.clone()]);
        let term1 = alpha1_xd.poly_add(&beta1_poly).unwrap();
        let h_term1 = c1_polynomial.poly_mul(&term1).unwrap();

        println!("H_term1(x): {:?}", h_term1);

        // C2(x) * (α2*x^(D-D2) + β2)
        let d2 = 3u32;
        let exponent = d - d2;
        let mut x_to_d_coeffs = vec![FieldElement::new(BigUint::zero(), prime.clone()).unwrap(); exponent as usize];
        x_to_d_coeffs.push(FieldElement::new(BigUint::one(), prime.clone()).unwrap());
        let x_to_d = Polynomial::new(x_to_d_coeffs);
        let alpha2_xd = x_to_d.poly_scale(&alpha2).unwrap();
        let beta2_poly = Polynomial::new(vec![beta2.clone()]);
        let term2 = alpha2_xd.poly_add(&beta2_poly).unwrap();
        let h_term2 = c2_polynomial.poly_mul(&term2).unwrap();

        println!("h_term2: {:?}", h_term2);

        // H(x) = C1(x) * (α1*x^(D-D1) + β1) + C2(x) * (α2*x^(D-D2) + β2)
        let h_polynomial = h_term1.poly_add(&h_term2).unwrap();

        println!("H(x): {:?}", h_polynomial);

        // Sampling outside the original domain
        // z = 8
        let z_value = FieldElement::new(BigUint::from_u32(8u32).unwrap(), prime.clone()).unwrap();

        // t(z)
        let t_z = t_polynomial.evaluate(&z_value).unwrap();
        println!("t(8) = {:?}", t_z); // t(8)

        // t(gz): （g = ω = 13）
        let g = omega.clone();
        let gz = g.mul(&z_value).unwrap();
        let t_gz = t_polynomial.evaluate(&gz).unwrap();
        println!("t(g * z) = t({}) = {:?}", gz.num, t_gz);

        // p1(z) = t(z) - 3
        let p1_z = t_z.sub(&three).unwrap();
        println!("p1(8) = t(8) - 3 = {:?}", p1_z);

        // Z1(z) = z - 1
        let z1_z = z_value.sub(&one).unwrap();
        println!("Z1(8) = 8 - 1 = {:?}", z1_z);

        // C1(z) = p1(z) / Z1(z)
        let z1_z_inv = z1_z.inv().unwrap();
        let c1_z = p1_z.mul(&z1_z_inv).unwrap();
        println!("C1(8) = p1(8) / Z1(8) = {:?}", c1_z);

        // C1(z) * (α1 * z^2 + β1)
        let z_squared = z_value.pow(2u32).unwrap();
        let term1_z = alpha1.mul(&z_squared).unwrap().add(&beta1).unwrap();
        let h_term1_z = c1_z.mul(&term1_z).unwrap();
        println!("H_term_1(8) = C1(8) * (α1 * 8^2 + β1) = {:?}", h_term1_z);

        // p2(z) = t(gz) - t(z)^2
        let t_z_squared = t_z.pow(2u32).unwrap();
        let p2_z = t_gz.sub(&t_z_squared).unwrap();
        println!("p2(8) = t(g * 8) - [t(8)]^2 = {:?}", p2_z);

        // z2(z) = z
        let z2_z = z_value.clone();
        println!("Z2(8) = {:?}", z2_z);

        // C2(z) = p2(z) / Z2(z)
        let z2_z_inv = z2_z.inv().unwrap();
        let c2_z = p2_z.mul(&z2_z_inv).unwrap();
        println!("C2(8) = p2(8) / Z2(8) = {:?}", c2_z);

        // C2(z) * (α2 * z + β2)
        let term2_z = alpha2.mul(&z_value).unwrap().add(&beta2).unwrap();
        let h_term2_z = c2_z.mul(&term2_z).unwrap();
        println!("H_term_2(8) = C2(8) * (α2 * 8 + β2) = {:?}", h_term2_z);

        // H(z)
        let h_z = h_term1_z.add(&h_term2_z).unwrap();
        println!("H(8) = H_term_1(8) + H_term_2(8) = {:?}", h_z);

        // H(8)
        let h_at_z = h_polynomial.evaluate(&z_value).unwrap();
        println!("H(8) (from H(x)) = {:?}", h_at_z);

        assert_eq!(h_z, h_at_z);

        // DEEP合成多項式

        // z = 8, gz = g * z, z^2 = z^2
        let z_value = FieldElement::new(BigUint::from_u32(8u32).unwrap(), prime.clone()).unwrap();
        let g = omega.clone();
        let gz_value = g.mul(&z_value).unwrap();
        let z_squared = z_value.pow(2u32).unwrap();

        // t(z), t(gz)
        let t_z = t_polynomial.evaluate(&z_value).unwrap();
        let t_gz = t_polynomial.evaluate(&gz_value).unwrap();

        // Term1: (t(x) - t(z)) / (x - z)
        let numerator1 = t_polynomial.poly_sub(&Polynomial::new(vec![t_z.clone()])).unwrap();
        let denominator1 = Polynomial::new(vec![
            z_value.clone().negate().unwrap(),
            FieldElement::new(BigUint::one(), prime.clone()).unwrap(),
        ]);

        let (term1_polynomial, remainder1) = numerator1.poly_div_rem(&denominator1).unwrap();

        println!("Term1 Polynomial: {:?}", term1_polynomial);
        println!("remainder1: {:?}", remainder1);

        // Term2: (t(x) - t(gz)) / (x - gz)
        let numerator2 = t_polynomial.poly_sub(&Polynomial::new(vec![t_gz.clone()])).unwrap();
        let denominator2 = Polynomial::new(vec![
            gz_value.clone().negate().unwrap(),
            FieldElement::new(BigUint::one(), prime.clone()).unwrap(),
        ]);

        let (term2_polynomial, remainder2) = numerator2.poly_div_rem(&denominator2).unwrap();

        println!("Term2 Polynomial: {:?}", term2_polynomial);
        println!("remainder2: {:?}", remainder2);

        // H1(x) と H2(x)
        let mut h1_coefficients = vec![];
        let mut h2_coefficients = vec![];

        // Coefficients of h_polynomial
        let h_coefficients = h_polynomial.coefficients.clone();

        for (i, coeff) in h_coefficients.iter().enumerate() {
            if i % 2 == 0 {
                // Even
                h1_coefficients.push(coeff.clone());
            } else {
                // Odd
                h2_coefficients.push(coeff.clone());
            }
        }

        let h1 = Polynomial::new(h1_coefficients);
        let h2 = Polynomial::new(h2_coefficients);

        println!("H1(x): {:?}", h1);
        println!("H2(x): {:?}", h2);

        // x^2
        let x_squared_poly = Polynomial::new(vec![
            FieldElement::new(BigUint::zero(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::zero(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::one(), prime.clone()).unwrap(),
        ]);

        // H1(x^2)
        let h1_x_squared = h1.compose(&x_squared_poly).unwrap();
        println!("H1(x^2): {:?}", h1_x_squared);

        // H2(x^2)
        let h2_x_squared = h2.compose(&x_squared_poly).unwrap();
        println!("H2(x^2): {:?}", h2_x_squared);

        // H1(z^2)
        let h1_z_squared = h1.evaluate(&z_squared).unwrap();
        println!("H1(z^2): {:?}", h1_z_squared);

        // H2(z^2)
        let h2_z_squared = h2.evaluate(&z_squared).unwrap();
        println!("H2(z^2): {:?}", h2_z_squared);

        // H1(x^2) - H1(z^2)
        let numerator3 = h1_x_squared.poly_sub(&Polynomial::new(vec![h1_z_squared.clone()])).unwrap();
        println!("Numerator3: {:?}", numerator3);

        // H2(x^2) - H2(z^2)
        let numerator4 = h2_x_squared.poly_sub(&Polynomial::new(vec![h2_z_squared.clone()])).unwrap();
        println!("Numerator4: {:?}", numerator4);

        // x - z
        let denominator = Polynomial::new(vec![
            z_value.clone().negate().unwrap(),
            FieldElement::new(BigUint::one(), prime.clone()).unwrap(),
        ]);

        println!("Denominator: {:?}", denominator);

        // Term3
        let (term3_polynomial, remainder3) = numerator3.poly_div_rem(&denominator).unwrap();

        println!("Term3 Polynomial: {:?}", term3_polynomial);
        println!("remainder3: {:?}", remainder3);

        // Term4
        let (term4_polynomial, remainder4) = numerator4.poly_div_rem(&denominator).unwrap();

        println!("Term4 Polynomial: {:?}", term4_polynomial);
        println!("remainder4: {:?}", remainder4);

        // γ_i
        let gamma1 = FieldElement::new(BigUint::one(), prime.clone()).unwrap();
        let gamma2 = FieldElement::new(BigUint::one(), prime.clone()).unwrap();
        let gamma3 = FieldElement::new(BigUint::one(), prime.clone()).unwrap();
        let gamma4 = FieldElement::new(BigUint::one(), prime.clone()).unwrap();

        // DEEP composition P0(x)
        let p0_polynomial = term1_polynomial.poly_scale(&gamma1).unwrap()
            .poly_add(&term2_polynomial.poly_scale(&gamma2).unwrap()).unwrap()
            .poly_add(&term3_polynomial.poly_scale(&gamma3).unwrap()).unwrap()
            .poly_add(&term4_polynomial.poly_scale(&gamma4).unwrap()).unwrap();

        println!("P0(x): {:?}", p0_polynomial);

        let expected_coefficients = vec![
            FieldElement::new(BigUint::from_u32(1u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(15u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(10u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(15u32).unwrap(), prime.clone()).unwrap(),
        ];

        let expected_p0 = Polynomial::new(expected_coefficients);

        println!("Expected P0(x): {:?}", expected_p0);

        assert_eq!(p0_polynomial, expected_p0);

    }
}
