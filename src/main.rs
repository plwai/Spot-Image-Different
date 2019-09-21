pub mod picture;

use std::env;
use std::time::Instant;

use crate::picture::Picture;

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut picture_1 = Picture::load(args.get(1).unwrap().to_string()).unwrap();
    let mut picture_2 = Picture::load(args.get(2).unwrap().to_string()).unwrap();
    let output_path = args.get(3).unwrap();

    let now = Instant::now();

    Picture::spot_different_custom_config(
        &mut picture_1,
        &mut picture_2,
        (5, 5),
        0.5,
        0.50,
        10.0,
        2,
        13,
        4,
        (255, 0, 0),
        output_path,
    );

    println!("Used {} seconds", now.elapsed().as_secs());
}
