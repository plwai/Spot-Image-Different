mod circle;
mod drawer;
mod point;
mod utils;

use image;
use image::*;

use std::collections::HashMap;

use self::drawer::draw_circle;
use self::point::Point;

pub struct Picture {
    img: DynamicImage,
    _img_path: String,
}

impl Picture {
    pub fn load(img_path: &str) -> Result<Self, std::io::Error> {
        match image::open(img_path) {
            Ok(img) => Ok(Picture {
                img,
                _img_path: img_path.to_string(),
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
                    draw_circle(
                        rhs_picture,
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
        eps: f32,
        min_points: usize,
        circle_offset: i32,
        circle_thickness: i32,
        circle_color: (u8, u8, u8),
        output_path: &str,
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
        let point_groups = utils::dbscan(&points, eps, min_points);
        let result = point_groups
            .iter()
            .fold(HashMap::new(), |mut acc, (point, group)| {
                let group = group.unwrap();
                let points = acc.entry(group).or_insert(vec![]);
                points.push(*point);
                acc
            });

        // Filter noise
        let result = result
            .into_iter()
            .filter(|(_group_id, points)| points.len() > min_points)
            .collect::<HashMap<i32, Vec<_>>>();

        // Find the smallest circle midpoint for every group of points
        let circles = circle::compute_smallest_circle(result);

        // Compute result
        for circle in circles {
            match circle {
                Some(c) => {
                    draw_circle(
                        rhs_picture,
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

        rhs_picture.img.save(output_path).unwrap();
    }
}
