use num_bigint::BigUint;
use num_traits::{Zero, FromPrimitive};
use crate::field::FieldElement;

#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    pub coefficients: Vec<FieldElement>,
}

// TODO: Use NTT, INTT
impl Polynomial {
    pub fn new(coefficients: Vec<FieldElement>) -> Self {
        Polynomial { coefficients }
    }

    // Addition between polynomials
    pub fn poly_add(&self, other: &Polynomial) -> Result<Polynomial, &'static str> {
        let len = usize::max(self.coefficients.len(), other.coefficients.len());
        let mut result_coeffs = Vec::with_capacity(len);
        let prime = self.coefficients[0].prime.clone();

        let zero = FieldElement::new(BigUint::zero(), prime.clone()).unwrap();
        for i in 0..len {
            let a = self.coefficients.get(i).unwrap_or(&zero);
            let b = other.coefficients.get(i).unwrap_or(&zero);
            let sum = a.add(b)?;
            result_coeffs.push(sum);
        }        

        Ok(Polynomial::new(result_coeffs))
    }

    // Subtraction between polynomials
    pub fn poly_sub(&self, other: &Polynomial) -> Result<Polynomial, &'static str> {
        let len = usize::max(self.coefficients.len(), other.coefficients.len());
        let mut result_coeffs = Vec::with_capacity(len);
        let prime = self.coefficients[0].prime.clone();

        let zero = FieldElement::new(BigUint::zero(), prime.clone()).unwrap();
        for i in 0..len {
            let a = self.coefficients.get(i).unwrap_or(&zero);
            let b = other.coefficients.get(i).unwrap_or(&zero);
            let diff = a.sub(b)?;
            result_coeffs.push(diff);
        }        

        Ok(Polynomial::new(result_coeffs))
    }

    // Multiplication between polynomials
    pub fn poly_mul(&self, other: &Polynomial) -> Result<Polynomial, &'static str> {
        let len = self.coefficients.len() + other.coefficients.len() - 1;
        let prime = self.coefficients[0].prime.clone();
        let zero = FieldElement::new(BigUint::zero(), prime.clone()).unwrap();
        let mut result_coeffs = vec![zero.clone(); len];

        for (i, coeff_a) in self.coefficients.iter().enumerate() {
            for (j, coeff_b) in other.coefficients.iter().enumerate() {
                let product = coeff_a.mul(coeff_b)?;
                let current = &result_coeffs[i + j];
                result_coeffs[i + j] = current.add(&product)?;
            }
        }

        Ok(Polynomial::new(result_coeffs))
    }

    // Scalar Multiplication of Polynomial
    pub fn poly_scale(&self, scalar: &FieldElement) -> Result<Polynomial, &'static str> {
        let scaled_coeffs = self.coefficients.iter().map(|c| c.mul(scalar)).collect::<Result<Vec<_>, _>>()?;
        Ok(Polynomial::new(scaled_coeffs))
    }

    // Division of polynomials and compute remainders
    pub fn poly_div_rem(&self, divisor: &Polynomial) -> Result<(Polynomial, Polynomial), &'static str> {
        let mut remainder = self.coefficients.clone();
        let mut quotient = Vec::new();
        let prime = self.coefficients[0].prime.clone();

        while remainder.len() >= divisor.coefficients.len() {
            let coeff = remainder.last().unwrap().div(&divisor.coefficients.last().unwrap())?;
            let degree = remainder.len() - divisor.coefficients.len();
            let zero = FieldElement::new(BigUint::zero(), prime.clone()).unwrap();

            let mut scaled_divisor = vec![zero.clone(); degree];
            scaled_divisor.extend(divisor.coefficients.iter().map(|c| c.mul(&coeff).unwrap()));

            remainder = remainder.iter()
                .zip(scaled_divisor.iter())
                .map(|(a, b)| a.sub(b).unwrap())
                .collect();

            while let Some(last) = remainder.last() {
                if last.num.is_zero() {
                    remainder.pop();
                } else {
                    break;
                }
            }

            quotient.insert(0, coeff);
        }

        Ok((Polynomial::new(quotient), Polynomial::new(remainder)))
    }

    // Evaluate Polynomial
    pub fn evaluate(&self, x: &FieldElement) -> Result<FieldElement, &'static str> {
        let result = evaluate_polynomial(&self.coefficients, x)?;
        Ok(result)
    }


    // Caluculate Composite polynomial
    pub fn compose(&self, other: &Polynomial) -> Result<Polynomial, String> {
        if self.coefficients.is_empty() {
            return Err("Cannot compose an empty polynomial.".to_string());
        }

        let prime = self.coefficients[0].prime.clone();
        let zero = FieldElement::new(BigUint::zero(), prime.clone()).unwrap();
        let mut result = Polynomial::new(vec![zero]);

        for coeff in self.coefficients.iter().rev() {
            result = result.poly_mul(other)?;
            result = result.poly_add(&Polynomial::new(vec![coeff.clone()]))?;
        }
        result.trim();
        Ok(result)
    }

    pub fn trim(&mut self) {
        while let Some(true) = self.coefficients.last().map(|c| c.is_zero()) {
            self.coefficients.pop();
        }
    }

    pub fn poly_pow(&self, exponent: u32, field_prime: &BigUint) -> Result<Self, &'static str> {
        let mut result = Polynomial::new(vec![FieldElement::one(field_prime)]);
        let mut base = self.clone();
        let mut exp = exponent;

        while exp > 0 {
            if exp % 2 == 1 {
                result = result.poly_mul(&base)?;
            }
            base = base.poly_mul(&base)?;
            exp /= 2;
        }

        Ok(result)
    }

    pub fn monomial(degree: usize, coefficient: FieldElement) -> Self {
        let mut coefficients = vec![FieldElement::zero(&coefficient.prime); degree];
        coefficients.push(coefficient);
        Polynomial::new(coefficients)
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|coeff| coeff.is_zero())
    }

    pub fn degree(&self) -> usize {
        let mut degree = self.coefficients.len();
        while degree > 0 && self.coefficients[degree - 1].is_zero() {
            degree -= 1;
        }
        degree.saturating_sub(1)
    }
}

// Evaluate polynomial.
pub fn evaluate_polynomial(coefficients: &[FieldElement], x: &FieldElement) -> Result<FieldElement, &'static str> {
    let mut result = FieldElement::new(BigUint::zero(), x.prime.clone()).unwrap();
    for (i, coeff) in coefficients.iter().enumerate() {
        let x_pow_i = x.pow(i as u32)?;
        let term = coeff.mul(&x_pow_i)?;
        result = result.add(&term)?;
    }
    Ok(result)
}

// Evaluate polynomial at multiple points.
pub fn evaluate_polynomial_at_points(coefficients: &[FieldElement], points: &[FieldElement]) -> Result<Vec<FieldElement>, &'static str> {
    let mut evaluations = Vec::new();
    for x in points {
        let eval = evaluate_polynomial(coefficients, x)?;
        evaluations.push(eval);
    }
    Ok(evaluations)
}

use std::convert::From;

impl From<Vec<FieldElement>> for Polynomial {
    fn from(coefficients: Vec<FieldElement>) -> Self {
        Polynomial::new(coefficients)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use num_traits::{FromPrimitive, Zero, One};
    use crate::field::FieldElement;
    use crate::ntt;

    #[test]
    fn test_intt_and_evaluation_2() {
        let prime = BigUint::from_u32(17u32).unwrap();

        let generator_value = BigUint::from_u32(3u32).unwrap();
        let generator_element = FieldElement::new(generator_value.clone(), prime.clone()).unwrap();
        let generator = generator_element.set_generator(prime.clone()).unwrap();

        let n = 4u32;
        let n_biguint = BigUint::from_u32(n).unwrap();

        let omega = generator.nth_generator(n_biguint.clone()).unwrap();
        let omega_inv = omega.inv().unwrap();

        let n_field = FieldElement::new(n_biguint.clone(), prime.clone()).unwrap();
        let n_inv = n_field.inv().unwrap();

        // [3, 9, 13, 16]
        let trace_sequence = vec![
            FieldElement::new(BigUint::from_u32(3u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(9u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap(),
        ];

        // INTT
        let coefficients = ntt::aggregate_intt(&trace_sequence, &omega_inv, &n_inv.clone()).unwrap();

        let expected_coefficients = vec![
            FieldElement::new(BigUint::from_u32(6u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(2u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap(),
        ];

        assert_eq!(coefficients, expected_coefficients);

        // Evaluation domain x_i = 3 * h^i
        // h = 9
        let h_value = BigUint::from_u32(9u32).unwrap();
        let h = FieldElement::new(h_value.clone(), prime.clone()).unwrap();

        // g = 3
        let g_value = BigUint::from_u32(3u32).unwrap();
        let g = FieldElement::new(g_value.clone(), prime.clone()).unwrap();

        // Shifted evaluation points
        let mut evaluation_points = Vec::new();
        for i in 0..8 {
            let h_i = h.pow(i as u32).unwrap();
            let x_i = g.mul(&h_i).unwrap();
            evaluation_points.push(x_i);
        }

        let expected_points = vec![
            FieldElement::new(BigUint::from_u32(3u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(10u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(5u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(11u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(14u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(7u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(12u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(6u32).unwrap(), prime.clone()).unwrap(),
        ];

        assert_eq!(evaluation_points, expected_points);

        let evaluations = evaluate_polynomial_at_points(&coefficients, &evaluation_points).unwrap();
        println!("評価結果: {:?}", evaluations);
        let expected = vec![
            FieldElement::new(BigUint::from_u32(15u32).unwrap(), prime.clone()).unwrap(), // P(3)
            FieldElement::new(BigUint::from_u32(4u32).unwrap(), prime.clone()).unwrap(),  // P(10)
            FieldElement::new(BigUint::from_u32(10u32).unwrap(), prime.clone()).unwrap(), // P(5)
            FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap(), // P(11)
            FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap(), // P(14)
            FieldElement::new(BigUint::from_u32(0u32).unwrap(), prime.clone()).unwrap(),  // P(7)
            FieldElement::new(BigUint::from_u32(0u32).unwrap(), prime.clone()).unwrap(),  // P(12)
            FieldElement::new(BigUint::from_u32(7u32).unwrap(), prime.clone()).unwrap(),  // P(6)
        ];

        assert_eq!(evaluations, expected);
    }

    #[test]
    fn test_intt_and_evaluation() {
        let prime = BigUint::from_u32(97u32).unwrap();

        let generator_value = BigUint::from_u32(5u32).unwrap();
        let generator_element = FieldElement::new(generator_value.clone(), prime.clone()).unwrap();
        let generator = generator_element.set_generator(prime.clone()).unwrap();

        let n = 4usize;
        let n_biguint = BigUint::from_usize(n).unwrap();

        let omega = generator.nth_generator(n_biguint.clone()).unwrap();

        println!("H2^2:{:?}", &omega.clone());
        let omega_inv = omega.inv().unwrap();

        let n_field = FieldElement::new(n_biguint.clone(), prime.clone()).unwrap();
        let n_inv = n_field.inv().unwrap();

        let trace_sequence = vec![
            FieldElement::new(BigUint::from_u32(3u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(9u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(81u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(62u32).unwrap(), prime.clone()).unwrap(),
        ];

        let coefficients = ntt::aggregate_intt(&trace_sequence, &omega_inv, &n_inv.clone()).unwrap();

        let expected_coefficients = vec![
            FieldElement::new(BigUint::from_u32(63u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(78u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(76u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(77u32).unwrap(), prime.clone()).unwrap(),
        ];

        assert_eq!(coefficients, expected_coefficients);

        // N = 32
        let n = 32usize;
        let n_biguint = BigUint::from_usize(n).unwrap();
        let h = generator.nth_generator(n_biguint.clone()).unwrap();
        let g_value = BigUint::from_u32(3u32).unwrap();
        let g = FieldElement::new(g_value.clone(), prime.clone()).unwrap();

        let mut evaluation_points = Vec::new();
        for i in 0..32 {
            let h_i = h.pow(i as u32).unwrap();
            let x_i = g.mul(&h_i).unwrap();
            evaluation_points.push(x_i);
        }

        println!("拡張かつシフト:{:?}", evaluation_points.clone());

        let evaluations = evaluate_polynomial_at_points(&coefficients, &evaluation_points).unwrap();
        println!("評価結果: {:?}", evaluations);
    }

    #[test]
    fn test_polynomial_compose() {
        let prime = BigUint::from_u32(17u32).unwrap();

        // f(x) = x + 1
        let one = FieldElement::new(BigUint::one(), prime.clone()).unwrap();
        let f = Polynomial::new(vec![one.clone(), one.clone()]);

        // g(x) = x^2
        let zero = FieldElement::new(BigUint::zero(), prime.clone()).unwrap();
        let g = Polynomial::new(vec![zero.clone(), zero.clone(), one.clone()]);

        // f(g(x))
        let composed = f.compose(&g).unwrap();
        println!("composed: {:?}", composed);

        // x^2 + 1
        let expected = Polynomial::new(vec![one.clone(), zero.clone(), one.clone()]);
        println!("expected: {:?}", expected);
    }

}