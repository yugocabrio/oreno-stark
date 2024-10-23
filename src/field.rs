use num_bigint::BigUint;
use num_traits::{Zero, one};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldElement {
    pub num: BigUint,
    pub prime: BigUint,
}

impl FieldElement {
    // Constructor to create a new FieldElement.
    pub fn new(num: BigUint, prime: BigUint) -> Result<Self, &'static str> {
        let num_mod_prime = num % &prime;
        Ok(FieldElement { num: num_mod_prime, prime })
    }

    // Ensure num and prime belong to the same F_p.
    pub fn check_same_prime(&self, other: &Self) -> Result<(), &'static str> {
        if self.prime != other.prime {
            return Err("Different Prime numbers");
        }
        Ok(())
    }

    // Add two FieldElements.
    pub fn add(&self, other: &Self) -> Result<Self, &'static str> {
        self.check_same_prime(other)?;
        let add_num = (&self.num + &other.num) % &self.prime;
        FieldElement::new(add_num, self.prime.clone())
    }

    // Subtract one FieldElement from another.
    pub fn sub(&self, other: &Self) -> Result<Self, &'static str> {
        self.check_same_prime(other)?;
        let sub_num = (&self.num + &self.prime - &other.num) % &self.prime;
        FieldElement::new(sub_num, self.prime.clone())
    }

    // Multiply two FieldElements.
    pub fn mul(&self, other: &Self) -> Result<Self, &'static str> {
        self.check_same_prime(other)?;
        let mul_num = (&self.num * &other.num) % &self.prime;
        FieldElement::new(mul_num, self.prime.clone())
    }

    // Calculate the multiplicative inverse of a FieldElement.
    pub fn inv(&self) -> Result<Self, &'static str> {
        if self.num.is_zero() {
            return Err("Cannot compute inverse of zero");
        }
        // Fermat's little theorem
        let extended = self.num.modpow(&(self.prime.clone() - BigUint::from(2u32)), &self.prime);
        FieldElement::new(extended, self.prime.clone())
    }

    // Divide one FieldElement by another.
    pub fn div(&self, other: &Self) -> Result<Self, &'static str> {
        let inv_other = other.inv()?;
        self.mul(&inv_other)
    }

    // Exponentiate a FieldElement by an integer exponent.
    pub fn pow(&self, exp: u32) -> Result<Self, &'static str> {
        let base = &self.num % &self.prime;
        let result_num = base.modpow(&BigUint::from(exp), &self.prime);
        Ok(FieldElement { num: result_num, prime: self.prime.clone() })
    }

    // Set the generator.
    pub fn set_generator(&self, prime: BigUint) -> Result<FieldElement, &'static str> {
        // TODO: Implement validation of generator as primitive root if needed.
        FieldElement::new(self.num.clone(), prime)
    }

    // nth_generator = generator^((P-1)/N)
    pub fn nth_generator(&self, nth: BigUint) -> Result<FieldElement, &'static str> {
        let exp = (&self.prime - 1u32) / nth;
        let nth_power = self.num.modpow(&exp, &self.prime);
        FieldElement::new(nth_power, self.prime.clone())
    }

    // Return the additive inverse of tthe FieldElement.
    pub fn negate(&self) -> Result<Self, &'static str> {
        let zero = FieldElement::new(BigUint::zero(), self.prime.clone()).unwrap();
        zero.sub(self)
    }

    // Check if the number is zero or not.
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    #[test]
    fn test_check_same_prime() {
        let prime1 = BigUint::from(7u32);
        let prime2 = BigUint::from(11u32);
        let a = FieldElement::new(BigUint::from(3u32), prime1.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(3u32), prime2.clone()).unwrap();
        assert!(a.check_same_prime(&b).is_err());
    }

    #[test]
    fn test_add() {
        let prime = BigUint::from(5u32);
        let a = FieldElement::new(BigUint::from(4u32), prime.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(1u32), prime.clone()).unwrap();
        assert_eq!(a.add(&b).unwrap().num, BigUint::from(0u32));
    }

    #[test]
    fn test_sub() {
        let prime = BigUint::from(5u32);
        let a = FieldElement::new(BigUint::from(2u32), prime.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        assert_eq!(a.sub(&b).unwrap().num, BigUint::from(4u32));
    }

    #[test]
    fn test_mul() {
        let prime = BigUint::from(5u32);
        let a = FieldElement::new(BigUint::from(2u32), prime.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        assert_eq!(a.mul(&b).unwrap().num, BigUint::from(1u32));
    }

    #[test]
    fn test_inv() {
        let prime = BigUint::from(5u32);
        let a = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        assert_eq!(a.inv().unwrap().num, BigUint::from(2u32));
    }

    #[test]
    fn test_div() {
        let prime = BigUint::from(5u32);
        let a = FieldElement::new(BigUint::from(4u32), prime.clone()).unwrap();
        let b = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        assert_eq!(a.div(&b).unwrap().num, BigUint::from(3u32));
    }

    #[test]
    fn test_negate() {
        let prime = BigUint::from(5u32);
        let a = FieldElement::new(BigUint::from(2u32), prime.clone()).unwrap();
        let expected_inverse = BigUint::from(3u32);
        assert_eq!(a.negate().unwrap().num, expected_inverse);
    }

    #[test]
    fn test_set_generator() {
        let prime = BigUint::from(17u32);
        let generator_value = FieldElement::new(BigUint::from(3u32), prime.clone()).unwrap();
        assert_eq!(generator_value.set_generator(prime).unwrap().num , BigUint::from(3u32));       
    }

    #[test]
    fn test_nth_generator_for_4th_root() {
        let prime = BigUint::from(17u32);
        let generator_value = BigUint::from(3u32);
        let generator = FieldElement::new(generator_value.clone(), prime.clone()).unwrap();

        let nth = BigUint::from(4u32);
        let nth_generator = generator.nth_generator(nth).unwrap();

        assert_eq!(nth_generator.num, BigUint::from(13u32));
    }

}
