mod picture;


use image;
use image::*;
use picture::Picture;

fn main() {
    let mut picture_1 = Picture::load("./left27.jpg".to_string()).unwrap();
    let mut picture_2 = Picture::load("./right27.jpg".to_string()).unwrap();

    //picture_1.spot_different(&mut picture_2, 0.4);
    picture_1.spot_different2(&mut picture_2);
}
