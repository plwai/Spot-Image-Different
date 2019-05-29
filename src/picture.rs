mod circle;
mod point;
pub mod utils;

use image;
use image::*;

use std::collections::HashMap;

pub use self::circle::Circle;
pub use self::point::Point;

pub struct Picture {
    img: DynamicImage,
    _img_path: String,
}

impl Picture {
    pub fn load(img_path: String) -> Result<Self, std::io::Error> {
        match image::open(&img_path) {
            Ok(img) => Ok(Picture {
                img,
                _img_path: img_path,
            }),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Path not found!",
            )),
        }
    }

    pub fn spot_different_default_config(&mut self, rhs_picture: &mut Picture) {
        // Get points which have high deviation
        let points =
            utils::get_diff_points(&mut self.img, &mut rhs_picture.img, (10, 10), 0.4, 0.9);

        // Cluster the points into groups
        let point_groups = utils::dbscan(&points, 20.0, 1);
        let result = point_groups
            .iter()
            .fold(HashMap::new(), |mut acc, (point, group)| {
                let group = group.unwrap();
                let points = acc.entry(group).or_insert(vec![]);
                points.push(*point);
                acc
            });

        // Find the smallest circle midpoint for every group of points
        let circles = circle::compute_smallest_circle(result);

        // Compute result
        for circle in circles {
            match circle {
                Some(c) => {
                    rhs_picture.draw_circle(
                        c.mid_point.x as i32,
                        c.mid_point.y as i32,
                        c.radius as i32 + 13,
                        c.radius as i32 + 17,
                        (255, 0, 0),
                    );
                }
                None => unimplemented!(),
            }
        }

        rhs_picture.img.save("answer.jpg").unwrap();
    }

    pub fn spot_different_custom_config(
        &mut self,
        rhs_picture: &mut Picture,
        block_dimension: (u32, u32),
        threshold: f32,
        threshold_delta_rate: f32,
        fps: f32,
        min_points: usize,
        circle_offset: i32,
        circle_thickness: i32,
        circle_color: (u8, u8, u8),
    ) {
        // Get points which have high deviation
        let points = utils::get_diff_points(
            &mut self.img,
            &mut rhs_picture.img,
            block_dimension,
            threshold,
            threshold_delta_rate,
        );

        // Cluster the points into groups
        let point_groups = utils::dbscan(&points, fps, min_points);
        let result = point_groups
            .iter()
            .fold(HashMap::new(), |mut acc, (point, group)| {
                let group = group.unwrap();
                let points = acc.entry(group).or_insert(vec![]);
                points.push(*point);
                acc
            });

        // Find the smallest circle midpoint for every group of points
        let circles = circle::compute_smallest_circle(result);

        // Compute result
        for circle in circles {
            match circle {
                Some(c) => {
                    rhs_picture.draw_circle(
                        c.mid_point.x as i32,
                        c.mid_point.y as i32,
                        c.radius as i32 + circle_offset,
                        c.radius as i32 + circle_offset + circle_thickness,
                        circle_color,
                    );
                }
                None => unimplemented!(),
            }
        }

        rhs_picture.img.save("answer.jpg").unwrap();
    }

    fn draw_pixel(&mut self, x: i32, y: i32, colour: (u8, u8, u8)) {
        let mut x = x;
        let mut y = y;

        if x >= self.img.dimensions().0 as i32 {
            x = self.img.dimensions().0 as i32 - 1;
        }

        if y >= self.img.dimensions().1 as i32 {
            y = self.img.dimensions().1 as i32 - 1;
        }

        if x < 0 {
            x = 0;
        }

        if y < 0 {
            y = 0;
        }

        self.img.put_pixel(
            x as u32,
            y as u32,
            image::Rgba([colour.0, colour.1, colour.2, 1]),
        );
    }

    // https://stackoverflow.com/questions/27755514/circle-with-thickness-drawing-algorithm
    fn draw_circle(
        &mut self,
        xc: i32,
        yc: i32,
        inner_radius: i32,
        outer_radius: i32,
        colour: (u8, u8, u8),
    ) {
        let mut xo = outer_radius;
        let mut xi = inner_radius;
        let mut y = 0;
        let mut erro = 1 - xo;
        let mut erri = 1 - xi;

        while xo >= y {
            self.x_line(xc + xi, xc + xo, yc + y, colour);
            self.y_line(xc + y, yc + xi, yc + xo, colour);
            self.x_line(xc - xo, xc - xi, yc + y, colour);
            self.y_line(xc - y, yc + xi, yc + xo, colour);
            self.x_line(xc - xo, xc - xi, yc - y, colour);
            self.y_line(xc - y, yc - xo, yc - xi, colour);
            self.x_line(xc + xi, xc + xo, yc - y, colour);
            self.y_line(xc + y, yc - xo, yc - xi, colour);

            y = y + 1;

            if erro < 0 {
                erro += 2 * y + 1;
            } else {
                xo = xo - 1;
                erro += 2 * (y - xo + 1);
            }

            if y > inner_radius {
                xi = y;
            } else {
                if erri < 0 {
                    erri += 2 * y + 1;
                } else {
                    xi = xi - 1;
                    erri += 2 * (y - xi + 1);
                }
            }
        }
    }

    fn x_line(&mut self, x1: i32, x2: i32, y: i32, colour: (u8, u8, u8)) {
        let mut x1 = x1;

        while x1 <= x2 {
            self.draw_pixel(x1, y, colour);
            x1 = x1 + 1;
        }
    }

    fn y_line(&mut self, x: i32, y1: i32, y2: i32, colour: (u8, u8, u8)) {
        let mut y1 = y1;

        while y1 <= y2 {
            self.draw_pixel(x, y1, colour);
            y1 = y1 + 1;
        }
    }
}