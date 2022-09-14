use rand::{CryptoRng, Rng, thread_rng};
use rand::distributions::{Distribution, Uniform, uniform::SampleUniform};
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;
use crate::{galois_fields, GaloisField};
use crate::galois_fields::gf4_number::GF4;

pub struct Context {
    rng: ThreadRng,
}

impl Context {
    pub fn new() -> Context {
        Context{
            rng: thread_rng(),
        }
    }

    pub fn random_vector<T: GaloisField>(&mut self, length: usize) -> Vec<T> {
        let mut v = Vec::new();
        for i in 0..length {
            v.push(T::generate_random(&mut self.rng));
        }
        v
    }

    pub fn random_error_vector<T: GaloisField>(&mut self, length: usize, weight: usize) -> Vec<T> {
        let mut v: Vec<T> = vec![T::generate_zero(); length];
        for i in 0..weight {
            let val: T = T::generate_random(&mut self.rng);

            v[i] = T::generate_random(&mut self.rng);
        }
        v.shuffle(&mut self.rng);
        v
    }
}

/*
impl<T: SampleUniform + ?Sized> FiniteFieldGenerator<T> {
    pub fn new(low_inclusive: T, high_inclusive: T) -> FiniteFieldGenerator<T> {
        FiniteFieldGenerator {
            distribution: Uniform::new_inclusive(low_inclusive, high_inclusive),
            rng: thread_rng(),
        }
    }

    pub fn get_random(&mut self) -> T {
        self.distribution.sample(&mut self.rng)
    }
}
*/