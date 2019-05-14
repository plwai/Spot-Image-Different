mod picture;

use std::fs::File;
use std::io::Write;

use image;
use image::*;

fn main() {
    let mut raw_data_file = File::create("./picture_raw.txt").unwrap();
    let img = image::open("./original.jpg").unwrap();

    for (x, y, data) in img.pixels() {
        let rgb = data.to_rgb();
        let [r, g, b] = rgb.data;

        write!(raw_data_file, "{} {} {} ", r, g, b).unwrap();

        if x == 0 && y != 0 {
            writeln!(raw_data_file, "").unwrap();
        }
    }
}
