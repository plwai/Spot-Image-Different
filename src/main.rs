use std::env;

use spot_difference::picture::Picture;

fn main() {
    let args: Vec<String> = env::args().collect();
    
    let mut picture_1 = Picture::load(args.get(1).to_string()).unwrap();
    let mut picture_2 = Picture::load(args.get(2).to_string()).unwrap();

    picture_1.spot_different_default_config(&mut picture_2);
}
