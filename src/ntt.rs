use num_bigint::BigUint;
use num_traits::Zero;
use crate::field::FieldElement;

pub fn ntt_recursive(a: &[FieldElement], omega: &FieldElement) -> Result<Vec<FieldElement>, &'static str> {
    let n = a.len();
    if n == 1 {
        return Ok(a.to_vec());
    }

    let omega_squared = omega.pow(2u32)?;
    
    let a_even: Vec<FieldElement> = a.iter().step_by(2).cloned().collect();
    let a_odd: Vec<FieldElement> = a.iter().skip(1).step_by(2).cloned().collect();

    let y_even = ntt_recursive(&a_even, &omega_squared)?;
    let y_odd = ntt_recursive(&a_odd, &omega_squared)?;

    let mut y = vec![FieldElement::new(BigUint::zero(), omega.prime.clone()).unwrap(); n];
    for k in 0..(n / 2) {
        let omega_k = omega.pow(k as u32)?;
        let t = y_odd[k].mul(&omega_k)?;
        y[k] = y_even[k].add(&t)?;
        y[k + n / 2] = y_even[k].sub(&t)?;
    }
    Ok(y)
}

pub fn aggregate_intt(a: &[FieldElement], omega_inv: &FieldElement, nth_inv: &FieldElement) -> Result<Vec<FieldElement>, &'static str> {
    let intt_result = intt_recursive(a, &omega_inv).unwrap();
    let scaled_result = scale(&intt_result, nth_inv).unwrap();
    Ok(scaled_result)
}

pub fn intt_recursive(a: &[FieldElement], omega_inv: &FieldElement) -> Result<Vec<FieldElement>, &'static str> {
    let n = a.len();
    if n == 1 {
        return Ok(a.to_vec());
    }

    let omega_inv_squared = omega_inv.pow(2u32)?;

    let a_even: Vec<FieldElement> = a.iter().step_by(2).cloned().collect();
    let a_odd: Vec<FieldElement> = a.iter().skip(1).step_by(2).cloned().collect();

    let y_even = intt_recursive(&a_even, &omega_inv_squared)?;
    let y_odd = intt_recursive(&a_odd, &omega_inv_squared)?;

    let mut y = vec![FieldElement::new(BigUint::zero(), omega_inv.prime.clone()).unwrap(); n];
    for k in 0..(n / 2) {
        let omega_inv_k = omega_inv.pow(k as u32)?;
        let t = y_odd[k].mul(&omega_inv_k)?;
        y[k] = y_even[k].add(&t)?;
        y[k + n / 2] = y_even[k].sub(&t)?;
    }
    Ok(y)
}

// Scaling in the INTT
pub fn scale(a: &[FieldElement], n_inv: &FieldElement) -> Result<Vec<FieldElement>, &'static str> {
    let scaled: Vec<FieldElement> = a.iter().map(|x| x.mul(n_inv)).collect::<Result<Vec<FieldElement>, &'static str>>()?;
    Ok(scaled)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use crate::field::FieldElement;

    fn setup() -> (BigUint, FieldElement, FieldElement, FieldElement, FieldElement) {
        // Set prime number
        let prime = BigUint::from_u32(17u32).unwrap();

        // The Generator
        let generator_value = BigUint::from_u32(3u32).unwrap();
        let generator_element = FieldElement::new(generator_value.clone(), prime.clone()).unwrap();
        let generator = generator_element.set_generator(prime.clone()).unwrap();

        // The number of the element in the trace as power of 2
        let nth = BigUint::from_u32(4u32).unwrap();

        // Nth generator and its inverse
        let nth_generator = generator.nth_generator(nth.clone()).unwrap();
        let nth_generator_inv = nth_generator.inv().unwrap();

        // Inverse of nth
        let nth_field = FieldElement::new(nth.clone(), prime.clone()).unwrap();
        let nth_inv = nth_field.inv().unwrap();

        (prime, nth_field, nth_generator, nth_generator_inv, nth_inv)
    }

    #[test]
    fn test_ntt() {
        let (prime, _nth_field, nth_generator, _nth_generator_inv, _nth_inv) = setup();

        // t(x) = 13x^3 + 2x^2 + 16x + 6
        let coefficients = vec![
            FieldElement::new(BigUint::from_u32(6u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(2u32).unwrap(), prime.clone()).unwrap(), 
            FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap(),
        ];

        let ntt_result = ntt_recursive(&coefficients, &nth_generator).unwrap();

        // Expected Evaluation form [3, 9, 13, 16]
        let expected = vec![
            FieldElement::new(BigUint::from_u32(3u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(9u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap(),
        ];

        assert_eq!(ntt_result, expected);
    }

    #[test]
    fn test_simple_intt() {
        let (prime, _nth_field, nth_generator, nth_generator_inv, nth_inv) = setup();

        // Evaluation form [3, 9, 13, 16]
        let evaluation_form = vec![
            FieldElement::new(BigUint::from_u32(3u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(9u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap(),
        ];

        // Expected Coefficients
        // t(x) = 13x^3 + 2x^2 + 16x + 6
        let expected_coefficients = vec![
            FieldElement::new(BigUint::from_u32(6u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(16u32).unwrap(), prime.clone()).unwrap(),
            FieldElement::new(BigUint::from_u32(2u32).unwrap(), prime.clone()).unwrap(), 
            FieldElement::new(BigUint::from_u32(13u32).unwrap(), prime.clone()).unwrap(),
        ];

        // let intt_result = intt_recursive(&evaluation_form, &nth_generator_inv).unwrap();
        // let scaled_result = scale(&intt_result, &nth_inv).unwrap();
        let scaled_result = aggregate_intt(&evaluation_form, &nth_generator_inv, &nth_inv);

        assert_eq!(scaled_result.unwrap(), expected_coefficients);
    }
}
