use crate::galois_fields::GaloisField;

mod galois_fields;
mod polynomials;
mod random;

fn main() {
    use galois_fields::gf4_number::GF4;
    use polynomials::polynomial::Polynomial;


    let p1 = Polynomial::new_from_coefficients(vec![
        GF4::Zero, GF4::One, GF4::Alpha, GF4::Alpha, GF4::AlphaPlusOne
    ]);
    let p2 = Polynomial::new_from_coefficients(vec![
        GF4::AlphaPlusOne, GF4::Zero, GF4::One
    ]);

    let p3 = p1.mul(&p2);
    println!("{:?}", p3);
}
