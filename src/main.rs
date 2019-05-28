mod picture;

use picture::Picture;

fn main() {
    let mut picture_1 = Picture::load("./left27.jpg".to_string()).unwrap();
    let mut picture_2 = Picture::load("./right27.jpg".to_string()).unwrap();

    picture_1.spot_different_default_config(&mut picture_2);
}
