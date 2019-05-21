mod picture;

use picture::Picture;

fn main() {
    let mut picture_1 = Picture::load("./left.jpg".to_string()).unwrap();
    let mut picture_2 = Picture::load("./right.jpg".to_string()).unwrap();

    picture_1.spot_different2(&mut picture_2);
}
