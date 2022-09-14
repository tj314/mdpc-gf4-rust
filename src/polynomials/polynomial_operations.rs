use crate::GaloisField;
use crate::polynomials::polynomial::Polynomial;

pub fn xgcd<T: GaloisField>(poly: &Polynomial<T>, modulus: &Polynomial<T>) -> (Option<Polynomial<T>>, Option<Polynomial<T>>) {
    if modulus.degree() <= poly.degree() || poly.is_zero() {
        (None, None)
    } else {
        let mut r_last = modulus.clone();
        let mut r_current = poly.clone();
        let mut t_last = Polynomial::<T>::new();
        let mut t_current = Polynomial::<T>::new_from_coefficients(vec![T::generate_zero()]);

        loop {
            if let Some((q_current, mod_current)) = r_last.div_mod(&r_current) {
                let t = t_last.sub(&q_current.mul(&t_current));
                r_last = r_current;
                r_current = mod_current;
                t_last = t_current;
                t_current = t;
                if r_current.is_one() {
                    return (Some(r_current), Some(t_current));
                } else if r_current.is_zero() {
                    return (Some(r_last), None);
                }
            } else {
                return (None, None);
            }
        }
    }
}

#[cfg(test)]
mod polynomial_operations_tests {
    use crate::galois_fields::gf4_number::GF4;
    use super::*;


    #[test]
    fn test_xgcd() {
        {
            let modulus = Polynomial::new_from_coefficients(vec![
                GF4::Zero, GF4::One, GF4::Zero, GF4::Alpha
            ]);
            let poly = Polynomial::new_from_coefficients(vec![
                GF4::One, GF4::Zero, GF4::One
            ]);

            let (maybe_gcd, maybe_inv) = xgcd(&poly, &modulus);
            assert!(maybe_gcd.is_some());
            assert!(maybe_inv.is_some());
            assert!(maybe_gcd.unwrap().is_one());
            let inv = maybe_inv.unwrap();
            assert_eq!(inv.degree(), 2);
            assert_eq!(inv.get_coefficient(0).unwrap(), GF4::One);
            assert_eq!(inv.get_coefficient(1).unwrap(), GF4::Zero);
            assert_eq!(inv.get_coefficient(2).unwrap(), GF4::AlphaPlusOne);
        }
        {
            let modulus = Polynomial::new_from_coefficients(vec![
                GF4::Zero, GF4::One, GF4::One, GF4::One
            ]);
            let poly = Polynomial::new_from_coefficients(vec![
                GF4::One, GF4::Alpha
            ]);

            let (maybe_gcd, maybe_inv) = xgcd(&poly, &modulus);
            assert!(maybe_gcd.is_some());
            assert!(maybe_inv.is_none());
            let gcd = maybe_gcd.unwrap();
            assert_eq!(gcd.degree(), 1);
            assert_eq!(gcd, poly);
        }
    }
}