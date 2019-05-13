use image;
use image::*;

use crate::algorithms::Algorithms;

struct Picture {
    img: DynamicImage,
    img_path: String,
}

impl Picture {
    fn load(img_path: String) -> Result<Self, std::io::Error> {
        match image::open(&img_path) {
            Ok(img) => Ok(Picture { img, img_path }),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Path not found!",
            )),
        }
    }

    fn spot_different(&mut self, rhs_picture: &mut Picture) {}
}

impl Algorithms for Picture {
    fn MRE() {}
    fn PSNR() {}
    fn SSIM() {}
}
