use rand::Rng;
use rand::seq::SliceRandom;
use crate::galois_fields::GaloisField;


#[derive(Eq, PartialEq, Debug, Clone)]
pub enum GF4 {
    Zero,
    One,
    Alpha,
    AlphaPlusOne,
}

static ADDITION: [[GF4; 4]; 4] = [
    [GF4::Zero, GF4::One, GF4::Alpha, GF4::AlphaPlusOne],
    [GF4::One, GF4::Zero, GF4::AlphaPlusOne, GF4::Alpha],
    [GF4::Alpha, GF4::AlphaPlusOne, GF4::Zero, GF4::One],
    [GF4::AlphaPlusOne, GF4::Alpha, GF4::One, GF4::Zero]
];

static MULTIPLICATION: [[GF4; 4]; 4] = [
    [GF4::Zero, GF4::Zero, GF4::Zero, GF4::Zero],
    [GF4::Zero, GF4::One, GF4::Alpha, GF4::AlphaPlusOne],
    [GF4::Zero, GF4::Alpha, GF4::AlphaPlusOne, GF4::One],
    [GF4::Zero, GF4::AlphaPlusOne, GF4::One, GF4::Alpha]
];

static DIVISION: [[GF4; 3]; 4] = [
    [GF4::Zero, GF4::Zero, GF4::Zero],
    [GF4::One, GF4::AlphaPlusOne, GF4::Alpha],
    [GF4::Alpha, GF4::One, GF4::AlphaPlusOne],
    [GF4::AlphaPlusOne, GF4::Alpha, GF4::One]
];

impl GF4 {
    pub fn to_number(&self) -> u8 {
        match self {
            GF4::Zero => 0,
            GF4::One => 1,
            GF4::Alpha => 2,
            GF4::AlphaPlusOne => 3,
        }
    }

    pub fn from_number(num: u8) -> Option<GF4> {
        match num {
            0 => Some(GF4::Zero),
            1 => Some(GF4::One),
            2 => Some(GF4::Alpha),
            3 => Some(GF4::AlphaPlusOne),
            _ => None,
        }
    }
}

impl GaloisField for GF4 {
    // type Output = GF4;
    
    fn generate_zero() -> GF4 {
        GF4::Zero
    }

    fn is_zero(&self) -> bool {
        self == &GF4::Zero
    }

    fn generate_one() -> GF4 {
        GF4::One
    }

    fn is_one(&self) -> bool {
        self == &GF4::One
    }

    fn generate_random<R: Rng + ?Sized>(rng: &mut R) -> GF4 {
        [GF4::Zero, GF4::One, GF4::Alpha, GF4::AlphaPlusOne].choose(rng).unwrap().clone()
    }

    fn add(&self, other: &GF4) -> GF4 {
        ADDITION[(self.to_number()) as usize][(other.to_number()) as usize].clone()
    }

    fn sub(&self, other: &GF4) -> GF4 {
        self.add(other)
    }

    fn mul(&self, other: &GF4) -> GF4 {
        MULTIPLICATION[(self.to_number()) as usize][(other.to_number()) as usize].clone()
    }

    fn div(&self, other: &GF4) -> Option<GF4> {
        if other.is_zero() {
            None
        } else {
            Some(DIVISION[(self.to_number()) as usize][(other.to_number()) as usize - 1].clone())
        }
    }
}

#[cfg(test)]
mod gf4_tests {
    use std::collections::HashMap;
    use super::*;

    #[test]
    fn test_gf4_to_number() {
        assert_eq!(GF4::Zero.to_number(), 0);
        assert_eq!(GF4::One.to_number(), 1);
        assert_eq!(GF4::Alpha.to_number(), 2);
        assert_eq!(GF4::AlphaPlusOne.to_number(), 3);
    }

    #[test]
    fn test_gf4_from_number() {
        assert_eq!(GF4::from_number(0).unwrap(), GF4::Zero);
        assert_eq!(GF4::from_number(1).unwrap(), GF4::One);
        assert_eq!(GF4::from_number(2).unwrap(), GF4::Alpha);
        assert_eq!(GF4::from_number(3).unwrap(), GF4::AlphaPlusOne);
        assert!(GF4::from_number(4).is_none());
    }

    #[test]
    fn test_gf4_generate_zero() {
        assert_eq!(GF4::generate_zero(), GF4::Zero);
    }

    #[test]
    fn test_gf4_is_zero() {
        assert!(GF4::Zero.is_zero());
        assert!(!GF4::One.is_zero());
        assert!(!GF4::Alpha.is_zero());
        assert!(!GF4::AlphaPlusOne.is_zero());
    }

    #[test]
    fn test_gf4_generate_one() {
        assert_eq!(GF4::generate_one(), GF4::One);
    }

    #[test]
    fn test_gf4_is_one() {
        assert!(!GF4::Zero.is_one());
        assert!(GF4::One.is_one());
        assert!(!GF4::Alpha.is_one());
        assert!(!GF4::AlphaPlusOne.is_one());
    }

    /*
    #[test]
    fn test_gf4_make_random_context() {
    }

    #[test]
    fn test_gf4_generate_random_uniform() {
        let mut ctx = FiniteFieldGenerator::new(0u8, 3u8);
        let generated = GF4::generate_random_uniform(&mut ctx);
        assert!(vec![GF4::Zero, GF4::One, GF4::Alpha, GF4::AlphaPlusOne].contains(&generated));

        // chi squared test of uniformity
        // there are 4 options => 3 degrees of freedom
        // probab alpha = 0.05
        // the value of chi squared statistic can therefore be precomputed
        // in this case, the value is 7.815
        // this value was sourced from https://medcalc.org/manual/chi-square-table.php
        // given we are testing uniformity and there are 4 options
        // we expect probability of generating a specific option to be 1/4
        // that means out of N samples, the expected count of each option is N*0.25

        let num_samples = 10000;
        let expected_per_option = num_samples as f64 / 4.0;
        let mut sample_counts = HashMap::<u8, usize>::new();

        sample_counts.insert(GF4::Zero.to_number(), 0);
        sample_counts.insert(GF4::One.to_number(), 0);
        sample_counts.insert(GF4::Alpha.to_number(), 0);
        sample_counts.insert(GF4::AlphaPlusOne.to_number(), 0);

        for _ in 0..num_samples {
            let generated = GF4::generate_random_uniform(&mut ctx).to_number();
            sample_counts.insert(generated, sample_counts.get(&generated).unwrap() + 1);
        }

        let mut stat = 0.0;
        for (key, value) in sample_counts {
            let diff_squared = (value as f64 - expected_per_option).powi(2);
            stat += diff_squared;
        }
        stat /= expected_per_option;

        assert!(stat < 7.815);
    }
    */

    #[test]
    fn test_gf4_add() {
        // 0 + x
        assert_eq!(GF4::Zero.add(&GF4::Zero), GF4::Zero);
        assert_eq!(GF4::Zero.add(&GF4::One), GF4::One);
        assert_eq!(GF4::Zero.add(&GF4::Alpha), GF4::Alpha);
        assert_eq!(GF4::Zero.add(&GF4::AlphaPlusOne), GF4::AlphaPlusOne);

        // 1 + x
        assert_eq!(GF4::One.add(&GF4::Zero), GF4::One);
        assert_eq!(GF4::One.add(&GF4::One), GF4::Zero);
        assert_eq!(GF4::One.add(&GF4::Alpha), GF4::AlphaPlusOne);
        assert_eq!(GF4::One.add(&GF4::AlphaPlusOne), GF4::Alpha);

        // a + x
        assert_eq!(GF4::Alpha.add(&GF4::Zero), GF4::Alpha);
        assert_eq!(GF4::Alpha.add(&GF4::One), GF4::AlphaPlusOne);
        assert_eq!(GF4::Alpha.add(&GF4::Alpha), GF4::Zero);
        assert_eq!(GF4::Alpha.add(&GF4::AlphaPlusOne), GF4::One);

        // (a+1) + x
        assert_eq!(GF4::AlphaPlusOne.add(&GF4::Zero), GF4::AlphaPlusOne);
        assert_eq!(GF4::AlphaPlusOne.add(&GF4::One), GF4::Alpha);
        assert_eq!(GF4::AlphaPlusOne.add(&GF4::Alpha), GF4::One);
        assert_eq!(GF4::AlphaPlusOne.add(&GF4::AlphaPlusOne), GF4::Zero);
    }

    #[test]
    fn test_gf4_sub() {
        // 0 - x
        assert_eq!(GF4::Zero.sub(&GF4::Zero), GF4::Zero);
        assert_eq!(GF4::Zero.sub(&GF4::One), GF4::One);
        assert_eq!(GF4::Zero.sub(&GF4::Alpha), GF4::Alpha);
        assert_eq!(GF4::Zero.sub(&GF4::AlphaPlusOne), GF4::AlphaPlusOne);

        // 1 - x
        assert_eq!(GF4::One.sub(&GF4::Zero), GF4::One);
        assert_eq!(GF4::One.sub(&GF4::One), GF4::Zero);
        assert_eq!(GF4::One.sub(&GF4::Alpha), GF4::AlphaPlusOne);
        assert_eq!(GF4::One.sub(&GF4::AlphaPlusOne), GF4::Alpha);

        // a - x
        assert_eq!(GF4::Alpha.sub(&GF4::Zero), GF4::Alpha);
        assert_eq!(GF4::Alpha.sub(&GF4::One), GF4::AlphaPlusOne);
        assert_eq!(GF4::Alpha.sub(&GF4::Alpha), GF4::Zero);
        assert_eq!(GF4::Alpha.sub(&GF4::AlphaPlusOne), GF4::One);

        // (a+1) - x
        assert_eq!(GF4::AlphaPlusOne.sub(&GF4::Zero), GF4::AlphaPlusOne);
        assert_eq!(GF4::AlphaPlusOne.sub(&GF4::One), GF4::Alpha);
        assert_eq!(GF4::AlphaPlusOne.sub(&GF4::Alpha), GF4::One);
        assert_eq!(GF4::AlphaPlusOne.sub(&GF4::AlphaPlusOne), GF4::Zero);
    }

    #[test]
    fn test_gf4_mul() {
        // 0 * x
        assert_eq!(GF4::Zero.mul(&GF4::Zero), GF4::Zero);
        assert_eq!(GF4::Zero.mul(&GF4::One), GF4::Zero);
        assert_eq!(GF4::Zero.mul(&GF4::Alpha), GF4::Zero);
        assert_eq!(GF4::Zero.mul(&GF4::AlphaPlusOne), GF4::Zero);

        // 1 * x
        assert_eq!(GF4::One.mul(&GF4::Zero), GF4::Zero);
        assert_eq!(GF4::One.mul(&GF4::One), GF4::One);
        assert_eq!(GF4::One.mul(&GF4::Alpha), GF4::Alpha);
        assert_eq!(GF4::One.mul(&GF4::AlphaPlusOne), GF4::AlphaPlusOne);

        // a * x
        assert_eq!(GF4::Alpha.mul(&GF4::Zero), GF4::Zero);
        assert_eq!(GF4::Alpha.mul(&GF4::One), GF4::Alpha);
        assert_eq!(GF4::Alpha.mul(&GF4::Alpha), GF4::AlphaPlusOne);
        assert_eq!(GF4::Alpha.mul(&GF4::AlphaPlusOne), GF4::One);

        // (a+1) * x
        assert_eq!(GF4::AlphaPlusOne.mul(&GF4::Zero), GF4::Zero);
        assert_eq!(GF4::AlphaPlusOne.mul(&GF4::One), GF4::AlphaPlusOne);
        assert_eq!(GF4::AlphaPlusOne.mul(&GF4::Alpha), GF4::One);
        assert_eq!(GF4::AlphaPlusOne.mul(&GF4::AlphaPlusOne), GF4::Alpha);
    }

    #[test]
    fn test_gf4_div() {
        // 0 / x
        assert!(GF4::Zero.div(&GF4::Zero).is_none());
        assert_eq!(GF4::Zero.div(&GF4::One).unwrap(), GF4::Zero);
        assert_eq!(GF4::Zero.div(&GF4::Alpha).unwrap(), GF4::Zero);
        assert_eq!(GF4::Zero.div(&GF4::AlphaPlusOne).unwrap(), GF4::Zero);

        // 1 / x
        assert!(GF4::One.div(&GF4::Zero).is_none());
        assert_eq!(GF4::One.div(&GF4::One).unwrap(), GF4::One);
        assert_eq!(GF4::One.div(&GF4::Alpha).unwrap(), GF4::AlphaPlusOne);
        assert_eq!(GF4::One.div(&GF4::AlphaPlusOne).unwrap(), GF4::Alpha);

        // a / x
        assert!(GF4::Alpha.div(&GF4::Zero).is_none());
        assert_eq!(GF4::Alpha.div(&GF4::One).unwrap(), GF4::Alpha);
        assert_eq!(GF4::Alpha.div(&GF4::Alpha).unwrap(), GF4::One);
        assert_eq!(GF4::Alpha.div(&GF4::AlphaPlusOne).unwrap(), GF4::AlphaPlusOne);

        // (a+1) / x
        assert!(GF4::AlphaPlusOne.div(&GF4::Zero).is_none());
        assert_eq!(GF4::AlphaPlusOne.div(&GF4::One).unwrap(), GF4::AlphaPlusOne);
        assert_eq!(GF4::AlphaPlusOne.div(&GF4::Alpha).unwrap(), GF4::Alpha);
        assert_eq!(GF4::AlphaPlusOne.div(&GF4::AlphaPlusOne).unwrap(), GF4::One);
    }
}