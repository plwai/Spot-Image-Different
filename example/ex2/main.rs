use std::time::Instant;

use spot_difference::picture::Picture;

fn main() {
    let mut picture_1 = Picture::load("./example/ex2/ex2_1.jpg").unwrap();
    let mut picture_2 = Picture::load("./example/ex2/ex2_2.jpg").unwrap();
    let output_path = "./example/ex2/Output.jpg";

    let now = Instant::now();

    Picture::spot_different_custom_config(
        &mut picture_1,
        &mut picture_2,
        (10, 10),
        0.4,
        0.9,
        20.0,
        1,
        13,
        4,
        (255, 0, 0),
        output_path,
    );

    println!("Used {} seconds", now.elapsed().as_secs());
}
