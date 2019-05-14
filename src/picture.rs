use image;
use image::*;

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

    fn ssim(img_1: DynamicImage, img_2 : DynamicImage) {
        let block_dimension = img_1.dimensions();
        let block_total_pixel = block_dimension.0 * block_dimension.1;

        let mut sum_x: u32 = 0;
        let mut sum_y: u32 = 0;

        for (x, y, data) in img_1.pixels() {
            let luma = data.to_luma();
            let [color_bit] = luma.data;

            sum_x = sum_x + color_bit as u32; 
        }

        for (x, y, data) in img_2.pixels() {
            let luma = data.to_luma();
            let [color_bit] = luma.data;

            sum_y = sum_y + color_bit as u32; 
        }

        let mu_x = sum_x / block_total_pixel;
        let mu_y = sum_y / block_total_pixel;

        //let (sigma_x, sigma_y, sigma_xy) = (0, 0, 0, 0);

        // SSIM Constants
        let bit : i32 = 2;
        let l = bit.pow(8) - 1;
        let k1 = 0.01;
        let k2 = 0.03;
    }
}
