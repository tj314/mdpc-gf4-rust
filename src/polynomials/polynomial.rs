use crate::galois_fields::GaloisField;
use crate::polynomials::polynomial_operations::xgcd;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial<T: GaloisField>{
    coefficients: Vec<T>,
}

fn remove_trailing_zeros<T: GaloisField>(v: &mut Vec<T>) {
    let mut count = 0usize;
    for item in v.iter().rev() {
        if item.is_zero() {
            count += 1;
        } else {
            break;
        }
    }
    let new_length = v.len().saturating_sub(count);
    v.truncate(new_length);

    if v.is_empty() {
        v.push(T::generate_zero());
    }
}


impl<T> Polynomial<T>
where T: GaloisField {

    fn shrink_to_degree(&mut self) {
        remove_trailing_zeros(&mut self.coefficients);
    }

    pub fn new() -> Polynomial<T> {
        Polynomial{
            coefficients: vec![T::generate_zero()]
        }
    }

    pub fn new_from_coefficients(coefficients: Vec<T>) -> Polynomial<T> {
        if coefficients.is_empty() {
            Polynomial::new()
        } else {
            let mut p = Polynomial{
                coefficients
            };
            p.shrink_to_degree();
            p
        }
    }

    pub fn is_zero(&self) -> bool {
        self.degree() == 0 && self.coefficients[0].is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.degree() == 0 && self.coefficients[0].is_one()
    }

    pub fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }

    pub fn get_coefficient(&self, i: usize) -> Option<T> {
        match self.coefficients.get(i) {
            Some(val) => Some(val.clone()),
            None => None
        }
    }

    pub fn add(&self, other: &Polynomial<T>) -> Polynomial<T> {
        let (longer, shorter) = if self.degree() >= other.degree() {
            (self, other)
        } else {
            (other, self)
        };

        let mut coefficients: Vec<T> = longer.coefficients.clone();
        for (i, item) in shorter.coefficients.iter().enumerate() {
            coefficients[i] = coefficients[i].add(item);
        }

        Polynomial::new_from_coefficients(coefficients)
    }

    pub fn sub(&self, other: &Polynomial<T>) -> Polynomial<T> {
        let (longer, shorter) = if self.degree() >= other.degree() {
            (self, other)
        } else {
            (other, self)
        };

        let mut coefficients: Vec<T> = longer.coefficients.clone();
        for (i, item) in shorter.coefficients.iter().enumerate() {
            coefficients[i] = coefficients[i].sub(item);
        }

        Polynomial::new_from_coefficients(coefficients)
    }

    pub fn mul(&self, other: &Polynomial<T>) -> Polynomial<T> {
        let mut coefficients = vec![T::generate_zero(); self.degree() + other.degree() + 1];
        for (i, item) in self.coefficients.iter().enumerate() {
            for (j, jtem) in other.coefficients.iter().enumerate() {
                // print!("{:?} + ({:?} * {:?}) = {:?} + {:?} = ", coefficients[i + j], item, jtem, coefficients[i + j], item.mul(jtem));
                coefficients[i + j] = coefficients[i + j].add(&item.mul(jtem));
                // println!("{:?}", coefficients[i + j]);
            }
            // println!();
        }
        Polynomial::new_from_coefficients(coefficients)
    }

    pub fn div_mod(&self, other: &Polynomial<T>) -> Option<(Polynomial<T>, Polynomial<T>)> {
        if self.degree() < other.degree() {
            Some((Polynomial::new(), self.clone()))
        } else if other.is_zero() {
            None
        } else {
            let mut current = self.coefficients.clone();
            let mut result = vec![T::generate_zero(); current.len()];
            while current.len() - 1 >= other.degree() {
                let d = current[current.len() - 1]
                    .div(&other.coefficients[other.coefficients.len() - 1])
                    .unwrap();
                let offset = current.len() - other.coefficients.len();
                result[offset] = d;
                for (i, item) in other.coefficients.iter().enumerate() {
                    current[i + offset] = current[i + offset].sub(&item.mul(&result[offset]));
                }
                remove_trailing_zeros(&mut current);
            }
            Some((Polynomial::new_from_coefficients(result), Polynomial::new_from_coefficients(current)))
        }
    }

    pub fn invert(&self, modulus: &Polynomial<T>) -> Option<Polynomial<T>> {
        let (_, maybe_inv) = xgcd(self, modulus);
        maybe_inv
    }
}

#[cfg(test)]
mod polynomial_tests {
    use crate::galois_fields::gf4_number::GF4;
    use super::*;


    /*
    Note: All tests in this module are performed using GF4, because GF4 has full test coverage,
    therefore GF4 is considered the reference implementation of GaloisField trait.
    */
    #[test]
    fn test_remove_trailing_zeros() {
        let mut input1 = vec![GF4::Zero, GF4::One, GF4::Alpha, GF4::Zero, GF4::Zero];
        let expected1 = vec![GF4::Zero, GF4::One, GF4::Alpha];

        let mut input2 = vec![GF4::Zero, GF4::One, GF4::Alpha, GF4::Zero, GF4::Zero, GF4::AlphaPlusOne];
        let expected2 = vec![GF4::Zero, GF4::One, GF4::Alpha, GF4::Zero, GF4::Zero, GF4::AlphaPlusOne];

        remove_trailing_zeros(&mut input1);
        remove_trailing_zeros(&mut input2);

        assert_eq!(input1, expected1);
        assert_eq!(input2, expected2);
    }

    #[test]
    fn test_polynomial_new() {
        let p: Polynomial<GF4> = Polynomial::new();
        assert_eq!(p.degree(), 0);
        assert_eq!(p.coefficients.len(), 1);
        assert_eq!(p.coefficients[0], GF4::Zero);
    }

    #[test]
    fn test_polynomial_new_from_coefficients() {
        let p1: Polynomial<GF4> = Polynomial::new_from_coefficients(vec![
            GF4::Zero, GF4::Zero, GF4::Zero
        ]);
        assert_eq!(p1.degree(), 0);
        assert_eq!(p1.coefficients.len(), 1);
        assert_eq!(p1.coefficients[0], GF4::Zero);

        let p1: Polynomial<GF4> = Polynomial::new_from_coefficients(vec![
            GF4::Zero, GF4::Zero, GF4::Zero, GF4::One, GF4::Zero
        ]);
        assert_eq!(p1.degree(), 3);
        assert_eq!(p1.coefficients.len(), 4);
        assert_eq!(p1.coefficients, vec![GF4::Zero, GF4::Zero, GF4::Zero, GF4::One]);
    }

    #[test]
    fn test_polynomial_is_zero() {
        let p1: Polynomial<GF4> = Polynomial::new();
        assert!(p1.is_zero());

        let p2: Polynomial<GF4> = Polynomial::new_from_coefficients(vec![GF4::One]);
        assert!(!p2.is_zero());
    }

    #[test]
    fn test_polynomial_is_one() {
        let p1: Polynomial<GF4> = Polynomial::new();
        assert!(!p1.is_one());

        let p2: Polynomial<GF4> = Polynomial::new_from_coefficients(vec![GF4::One]);
        assert!(p2.is_one());
    }

    #[test]
    fn test_polynomial_degree() {
        let p1: Polynomial<GF4> = Polynomial::new();
        let p2 = Polynomial::new_from_coefficients(vec![GF4::One, GF4::Zero, GF4::AlphaPlusOne, GF4::Zero]);
        let p3 = Polynomial::new_from_coefficients(vec![GF4::One, GF4::Zero, GF4::AlphaPlusOne]);
        assert_eq!(p1.degree(), 0);
        assert_eq!(p2.degree(), 2);
        assert_eq!(p3.degree(), 2);
    }

    #[test]
    fn test_polynomial_get_coefficient() {
        let p = Polynomial::new_from_coefficients(
            vec![GF4::Zero, GF4::One, GF4::Alpha, GF4::AlphaPlusOne, GF4::Zero]
        );
        assert!(p.get_coefficient(0).is_some());
        assert_eq!(p.get_coefficient(0).unwrap(), GF4::Zero);

        assert!(p.get_coefficient(1).is_some());
        assert_eq!(p.get_coefficient(1).unwrap(), GF4::One);

        assert!(p.get_coefficient(2).is_some());
        assert_eq!(p.get_coefficient(2).unwrap(), GF4::Alpha);

        assert!(p.get_coefficient(3).is_some());
        assert_eq!(p.get_coefficient(3).unwrap(), GF4::AlphaPlusOne);

        assert!(p.get_coefficient(4).is_none());
    }

    #[test]
    fn test_polynomial_add() {
        let p1 = Polynomial::new_from_coefficients(vec![
            GF4::Zero, GF4::One, GF4::Alpha, GF4::AlphaPlusOne
        ]);
        let p2 = Polynomial::new_from_coefficients(vec![
            GF4::One, GF4::One
        ]);
        let p3: Polynomial<GF4> = Polynomial::new();

        let p4 = p1.add(&p2);
        let p5 = p1.add(&p3);
        let p6 = p2.add(&p1);

        assert_eq!(p4.coefficients, vec![GF4::One, GF4::Zero, GF4::Alpha, GF4::AlphaPlusOne]);
        assert_eq!(p5.coefficients, p1.coefficients);
        assert_eq!(p6.coefficients, vec![GF4::One, GF4::Zero, GF4::Alpha, GF4::AlphaPlusOne]);
    }

    #[test]
    fn test_polynomial_sub() {
        let p1 = Polynomial::new_from_coefficients(vec![
            GF4::Zero, GF4::One, GF4::Alpha, GF4::AlphaPlusOne
        ]);
        let p2 = Polynomial::new_from_coefficients(vec![
            GF4::One, GF4::One
        ]);
        let p3: Polynomial<GF4> = Polynomial::new();

        let p4 = p1.sub(&p2);
        let p5 = p1.sub(&p3);
        let p6 = p2.sub(&p1);

        assert_eq!(p4.coefficients, vec![GF4::One, GF4::Zero, GF4::Alpha, GF4::AlphaPlusOne]);
        assert_eq!(p5.coefficients, p1.coefficients);
        assert_eq!(p6.coefficients, vec![GF4::One, GF4::Zero, GF4::Alpha, GF4::AlphaPlusOne]);
    }

    #[test]
    fn test_polynomial_mul() {
        let p1 = Polynomial::new_from_coefficients(vec![
            GF4::Zero, GF4::One, GF4::Alpha, GF4::Alpha, GF4::AlphaPlusOne
        ]);
        let p2 = Polynomial::new_from_coefficients(vec![
            GF4::AlphaPlusOne, GF4::Zero, GF4::One
        ]);

        let p3 = p1.mul(&p2);
        assert_eq!(p3.coefficients, vec![
            GF4::Zero, GF4::AlphaPlusOne, GF4::One, GF4::Zero, GF4::Zero, GF4::Alpha, GF4::AlphaPlusOne
        ]);
    }

    #[test]
    fn test_polynomial_div_mod() {
        let p1 = Polynomial::<GF4>::new();
        let p2 = Polynomial::new_from_coefficients(vec![
            GF4::One, GF4::Zero, GF4::One, GF4::Alpha, GF4::AlphaPlusOne
        ]);
        let p3 = Polynomial::new_from_coefficients(vec![GF4::One, GF4::Alpha]);

        let p4 = p2.div_mod(&p1);
        assert!(p4.is_none());

        let p5 = p2.div_mod(&p3);
        assert!(p5.is_some());
        let (p5_div, p5_mod) = p5.unwrap();
        assert_eq!(p5_div.coefficients, vec![GF4::Alpha, GF4::AlphaPlusOne, GF4::Zero, GF4::Alpha]);
        assert_eq!(p5_mod.coefficients, vec![GF4::AlphaPlusOne]);

        let p6 = p3.div_mod(&p2);
        assert!(p6.is_some());
        let (p6_div, p6_mod) = p6.unwrap();
        assert!(p6_div.is_zero());
        assert_eq!(p6_mod.coefficients, p3.coefficients);
    }

    #[test]
    fn test_polynomial_invert() {
        {
            let p = Polynomial::new_from_coefficients(vec![
                GF4::One, GF4::Zero, GF4::One
            ]);
            let m = Polynomial::new_from_coefficients(vec![
                GF4::Zero, GF4::One, GF4::Zero, GF4::Alpha
            ]);
            let inv = p.invert(&m);
            assert!(inv.is_some());
            assert_eq!(inv.unwrap().coefficients, vec![GF4::One, GF4::Zero, GF4::AlphaPlusOne]);
        }
        {
            let p = Polynomial::new_from_coefficients(vec![
                GF4::One, GF4::Alpha
            ]);
            let m = Polynomial::new_from_coefficients(vec![
                GF4::Zero, GF4::One, GF4::One, GF4::One
            ]);
            let inv = p.invert(&m);
            assert!(inv.is_none());
        }
    }
}