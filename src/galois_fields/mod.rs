use std::fmt::Debug;
use rand::distributions::{Distribution, Uniform};
use rand::Rng;

pub mod gf4_number;


pub trait GaloisField: Clone + Eq + PartialEq<Self> + Debug {
    fn generate_zero() -> Self;
    fn is_zero(&self) -> bool;
    fn generate_one() -> Self;
    fn is_one(&self) -> bool;
    fn generate_random<R: Rng + ?Sized>(rng: &mut R) -> Self;
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn div(&self, other: &Self) -> Option<Self>;
}