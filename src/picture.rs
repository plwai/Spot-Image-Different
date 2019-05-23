use image;
use image::*;

use std::collections::HashMap;

pub struct Picture {
    img: DynamicImage,
    img_path: String,
}

impl Picture {
    pub fn load(img_path: String) -> Result<Self, std::io::Error> {
        match image::open(&img_path) {
            Ok(img) => Ok(Picture { img, img_path }),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Path not found!",
            )),
        }
    }

    // SSIM comparison with only one block loop
    pub fn spot_different(&mut self, rhs_picture: &mut Picture, threshold: f32) {
        let block_dimension = (10 as u32, 10 as u32);

        let col_num = (self.img.dimensions().0 as f32 / block_dimension.0 as f32).ceil() as u32;
        let row_num = (self.img.dimensions().1 as f32 / block_dimension.1 as f32).ceil() as u32;

        let mut img_block_map = HashMap::new();

        for row in 0..row_num {
            for col in 0..col_num {
                let one_dimension = row * col_num + col;
                let current_origin = (block_dimension.0 * col, block_dimension.1 * row);

                let block_img_1 = self.img.crop(
                    current_origin.0,
                    current_origin.1,
                    block_dimension.0,
                    block_dimension.1,
                );
                let block_img_2 = rhs_picture.img.crop(
                    current_origin.0,
                    current_origin.1,
                    block_dimension.0,
                    block_dimension.1,
                );

                let ssim = Picture::ssim_rgb(&block_img_1, &block_img_2);

                if ssim < threshold {
                    println!("{} {} {}", current_origin.0, current_origin.1, ssim);
                    img_block_map.insert(one_dimension, true);
                } else {
                    img_block_map.insert(one_dimension, false);
                }
            }
        }

        let mut skip = vec![];

        for x in img_block_map.iter() {
            if skip.contains(x.0) {
                continue;
            }

            if *x.1 == false {
                continue;
            }

            let row = x.0 / col_num;
            let col = x.0 % col_num;

            let mut current_origin = (block_dimension.0 * col, block_dimension.1 * row);

            // Check right and bottom to merge block
            if col + 1 < col_num {
                let right = img_block_map.get(&(row * col_num + col + 1)).unwrap();

                if *right == true {
                    current_origin.0 = current_origin.0 + block_dimension.0;
                    skip.push(row * col_num + col + 1);
                }
            }

            if row + 1 < row_num {
                let btm = img_block_map.get(&((row + 1) * col_num + col)).unwrap();

                if *btm == true {
                    current_origin.1 = current_origin.1 + block_dimension.1;
                    skip.push((row + 1) * col_num + col);
                }
            }

            let crop_origin_x = match current_origin.0.overflowing_sub(20) {
                (x, false) => x,
                _ => 0,
            };

            let crop_origin_y = match current_origin.1.overflowing_sub(20) {
                (x, false) => x,
                _ => 0,
            };

            let block_img_1 = self.img.crop(crop_origin_x, crop_origin_y, 40, 40);
            let block_img_2 = rhs_picture.img.crop(crop_origin_x, crop_origin_y, 40, 40);
            block_img_1
                .save(format!(
                    "save_1_{}_{}_L.jpg",
                    current_origin.0, current_origin.1
                ))
                .unwrap();
            block_img_2
                .save(format!(
                    "save_2_{}_{}_L.jpg",
                    current_origin.0, current_origin.1
                ))
                .unwrap();
        }
    }

    pub fn spot_different2(&mut self, rhs_picture: &mut Picture) {
        let result =
            Picture::get_diff_points(&mut self.img, &mut rhs_picture.img, (10, 10), 0.4, 0.9);

        println!("");

        for p in result.iter() {
            print!("{} {}, ", p.x, p.y);
        }

        let result = Picture::dbscan(&result, 2.0, 1);
        let result = result
            .iter()
            .fold(HashMap::new(), |mut acc, (point, group)| {
                let group = group.unwrap();
                let points = acc.entry(group).or_insert(vec![]);
                points.push(point);
                acc
            });

        for (group, points) in result.iter() {
            println!("");
            println!("");
            println!("Group {}:", group);
            for p in points.iter() {
                print!("{} {}, ", p.x, p.y);
            }
        }
    }

    // Compare using SSIM algorithm from large blocks to small blocks
    // threshold < 1.0
    // threshold_delta_rate between 0 - 1.0
    fn get_diff_points(
        img_1: &mut DynamicImage,
        img_2: &mut DynamicImage,
        block_dimension: (u32, u32),
        threshold: f32,
        threshold_delta_rate: f32,
    ) -> Vec<Point> {
        let col_num = (img_1.dimensions().0 as f32 / block_dimension.0 as f32).ceil() as u32;
        let row_num = (img_1.dimensions().1 as f32 / block_dimension.1 as f32).ceil() as u32;

        let mut local_points = vec![];

        for row in 0..row_num {
            for col in 0..col_num {
                let current_origin = (block_dimension.0 * col, block_dimension.1 * row);

                let mut block_img_1 = img_1.crop(
                    current_origin.0,
                    current_origin.1,
                    block_dimension.0,
                    block_dimension.1,
                );
                let mut block_img_2 = img_2.crop(
                    current_origin.0,
                    current_origin.1,
                    block_dimension.0,
                    block_dimension.1,
                );

                let ssim = Picture::ssim_rgb(&block_img_1, &block_img_2);

                if ssim < threshold {
                    println!("{} {} {}", current_origin.0, current_origin.1, ssim);

                    if block_dimension != (1, 1) {
                        let inner_block_dimension = (
                            (block_dimension.0 as f32 / 2.0).floor() as u32,
                            (block_dimension.1 as f32 / 2.0).floor() as u32,
                        );

                        let inner_points = Self::get_diff_points(
                            &mut block_img_1,
                            &mut block_img_2,
                            inner_block_dimension,
                            threshold * threshold_delta_rate,
                            threshold_delta_rate,
                        );

                        for point in inner_points.iter() {
                            local_points.push(Point::new(
                                current_origin.0 + point.x,
                                current_origin.1 + point.y,
                            ));
                        }
                    } else {
                        local_points.push(Point::new(current_origin.0, current_origin.1));
                    }
                }
            }
        }

        local_points
    }

    fn dbscan(
        points: &Vec<Point>,
        eps: f32,
        minPts: usize,
    ) -> HashMap<Point, std::option::Option<i32>> {
        let mut cluster_id = 0;
        let mut points_map = points.iter().fold(HashMap::new(), |mut acc, p| {
            acc.insert(*p, None);
            acc
        });

        for point in points.iter() {
            if points_map[point].is_some() {
                continue;
            }

            let mut neighbors = points
                .iter()
                .filter(|p| Self::calculate_distance(&p, &point) <= eps)
                .collect::<Vec<_>>();

            if neighbors.len() < minPts {
                // -1 as noise
                points_map.insert(*point, Some(-1));
            }

            cluster_id = cluster_id + 1;
            points_map.insert(*point, Some(cluster_id));

            let mut start = 0;

            loop {
                for index in start..neighbors.len() {
                    let sp = neighbors[index];

                    if points_map[sp] == Some(-1) {
                        points_map.insert(*sp, Some(cluster_id));
                    }

                    if points_map[sp].is_some() {
                        continue;
                    }

                    points_map.insert(*sp, Some(cluster_id));

                    let seed_neighbors = points
                        .iter()
                        .filter(|p| {
                            Self::calculate_distance(p, sp) <= eps && !neighbors.contains(p)
                        })
                        .collect::<Vec<_>>();

                    if seed_neighbors.len() >= minPts {
                        neighbors.extend(seed_neighbors);
                    }

                    start = index;
                }

                if start == neighbors.len() - 1 {
                    break;
                }

                /*let neighbors = neighbors.iter().fold(vec!(), |acc, sp| {
                    if points_map[sp] == Some(-1) {
                        points_map.insert(**sp, Some(cluster_id));
                    }

                    if points_map[sp].is_some() {
                        acc.push(*sp);
                        return acc
                    }

                    points_map.insert(**sp, Some(cluster_id));

                    let mut seed_neighbors = points.iter().filter(|p| {
                        Self::calculate_distance(p, sp) <= eps && !neighbors.contains(p)
                    }).collect::<Vec<_>>();

                    if seed_neighbors.len() >= minPts {
                        acc.append(&mut seed_neighbors);
                    }

                    acc
                });*/
            }
        }

        points_map
    }

    fn calculate_distance(point_1: &Point, point_2: &Point) -> f32 {
        (((point_1.x as i32 - point_2.x as i32).abs().pow(2)
            + (point_1.y as i32 - point_2.y as i32).abs().pow(2)) as f32)
            .sqrt()
    }


    // SSIM using Luma data
    fn ssim(img_1: &DynamicImage, img_2: &DynamicImage) -> f32 {
        let img_1 = img_1.grayscale();
        let img_2 = img_2.grayscale();

        let block_dimension = img_1.dimensions();
        let block_total_pixel = (block_dimension.0 * block_dimension.1) as i64;

        let sum_x = img_1
            .pixels()
            .map(|(_x, _y, data)| {
                let luma = data.to_luma();
                let [color_bit] = luma.data;
                let color_bit = color_bit as i64;

                color_bit
            })
            .sum::<i64>();

        let sum_y = img_2
            .pixels()
            .map(|(_x, _y, data)| {
                let luma = data.to_luma();
                let [color_bit] = luma.data;
                let color_bit = color_bit as i64;

                color_bit
            })
            .sum::<i64>();

        let mu_x = sum_x / block_total_pixel;
        let mu_y = sum_y / block_total_pixel;

        let mut total_diff_mu_xy = 0;

        for x in 0..block_dimension.0 {
            for y in 0..block_dimension.1 {
                let [img_1_data] = img_1.get_pixel(x, y).to_luma().data;
                let [img_2_data] = img_2.get_pixel(x, y).to_luma().data;

                total_diff_mu_xy =
                    total_diff_mu_xy + ((img_1_data as i64 - mu_x) * (img_2_data as i64 - mu_y));
            }
        }


        let sigma_x = img_1
            .pixels()
            .map(|(_x, _y, data)| {
                let luma = data.to_luma();
                let [color_bit] = luma.data;
                let color_bit = color_bit as i64;

                (mu_x - color_bit).pow(2)
            })
            .sum::<i64>()
            / block_total_pixel;

        let sigma_y = img_2
            .pixels()
            .map(|(_x, _y, data)| {
                let luma = data.to_luma();
                let [color_bit] = luma.data;
                let color_bit = color_bit as i64;

                (mu_y - color_bit).pow(2)
            })
            .sum::<i64>()
            / block_total_pixel;

        let sigma_xy = total_diff_mu_xy / (block_total_pixel - 1);

        // SSIM Constants
        let bit: i32 = 2;
        let l = bit.pow(8) - 1;
        let k1 = 0.01;
        let k2 = 0.03;
        let c1 = (k1 * l as f32).powf(2.0);
        let c2 = (k2 * l as f32).powf(2.0);

        (((2 * mu_x * mu_y) as f32 + c1) * ((2 * sigma_xy) as f32 + c2))
            / (((mu_x.pow(2) + mu_y.pow(2)) as f32 + c1) * ((sigma_x + sigma_y) as f32 + c2))
    }

    // SSIM using RGB data
    fn ssim_rgb(img_1: &DynamicImage, img_2: &DynamicImage) -> f32 {
        let block_dimension = img_1.dimensions();
        let block_total_pixel = (block_dimension.0 * block_dimension.1) as i64;

        let sum_x = img_1
            .pixels()
            .map(|(_x, _y, data)| {
                let rgb = data.to_rgb();
                let [r, g, b] = rgb.data;
                let r = r as i64;
                let g = g as i64;
                let b = b as i64;

                (r, g, b)
            })
            .fold((0, 0, 0), |mut acc, x| {
                acc.0 = acc.0 + x.0;
                acc.1 = acc.1 + x.1;
                acc.2 = acc.2 + x.2;

                acc
            });

        let sum_y = img_2
            .pixels()
            .map(|(_x, _y, data)| {
                let rgb = data.to_rgb();
                let [r, g, b] = rgb.data;
                let r = r as i64;
                let g = g as i64;
                let b = b as i64;

                (r, g, b)
            })
            .fold((0, 0, 0), |mut acc, x| {
                acc.0 = acc.0 + x.0;
                acc.1 = acc.1 + x.1;
                acc.2 = acc.2 + x.2;

                acc
            });

        let mu_x = (
            sum_x.0 / block_total_pixel,
            sum_x.1 / block_total_pixel,
            sum_x.2 / block_total_pixel,
        );
        let mu_y = (
            sum_y.0 / block_total_pixel,
            sum_y.1 / block_total_pixel,
            sum_y.2 / block_total_pixel,
        );

        let mut total_diff_mu_xy = (0, 0, 0);

        for x in 0..block_dimension.0 {
            for y in 0..block_dimension.1 {
                let [img_1_r, img_1_g, img_1_b] = img_1.get_pixel(x, y).to_rgb().data;
                let [img_2_r, img_2_g, img_2_b] = img_2.get_pixel(x, y).to_rgb().data;

                total_diff_mu_xy.0 =
                    total_diff_mu_xy.0 + ((img_1_r as i64 - mu_x.0) * (img_2_r as i64 - mu_y.0));
                total_diff_mu_xy.1 =
                    total_diff_mu_xy.1 + ((img_1_g as i64 - mu_x.1) * (img_2_g as i64 - mu_y.1));
                total_diff_mu_xy.2 =
                    total_diff_mu_xy.2 + ((img_1_b as i64 - mu_x.2) * (img_2_b as i64 - mu_y.2));
            }
        }


        let sigma_x_mul_n = img_1
            .pixels()
            .map(|(_x, _y, data)| {
                let rgb = data.to_rgb();
                let [r, g, b] = rgb.data;
                let r = r as i64;
                let g = g as i64;
                let b = b as i64;

                (
                    (mu_x.0 - r).pow(2),
                    (mu_x.1 - g).pow(2),
                    (mu_x.2 - b).pow(2),
                )
            })
            .fold((0, 0, 0), |mut acc, x| {
                acc.0 = acc.0 + x.0;
                acc.1 = acc.1 + x.1;
                acc.2 = acc.2 + x.2;

                acc
            });

        let sigma_y_mul_n = img_2
            .pixels()
            .map(|(_x, _y, data)| {
                let rgb = data.to_rgb();
                let [r, g, b] = rgb.data;
                let r = r as i64;
                let g = g as i64;
                let b = b as i64;

                (
                    (mu_x.0 - r).pow(2),
                    (mu_x.1 - g).pow(2),
                    (mu_x.2 - b).pow(2),
                )
            })
            .fold((0, 0, 0), |mut acc, x| {
                acc.0 = acc.0 + x.0;
                acc.1 = acc.1 + x.1;
                acc.2 = acc.2 + x.2;

                acc
            });

        let sigma_x = (
            sigma_x_mul_n.0 / block_total_pixel,
            sigma_x_mul_n.1 / block_total_pixel,
            sigma_x_mul_n.2 / block_total_pixel,
        );

        let sigma_y = (
            sigma_y_mul_n.0 / block_total_pixel,
            sigma_y_mul_n.1 / block_total_pixel,
            sigma_y_mul_n.2 / block_total_pixel,
        );

        let sigma_xy = (
            total_diff_mu_xy.0 / (block_total_pixel),
            total_diff_mu_xy.1 / (block_total_pixel),
            total_diff_mu_xy.2 / (block_total_pixel),
        );

        // SSIM Constants
        let bit: i32 = 2;
        let l = bit.pow(8) - 1;
        let k1 = 0.01;
        let k2 = 0.03;
        let c1 = (k1 * l as f32).powf(2.0);
        let c2 = (k2 * l as f32).powf(2.0);

        let ssim_r = (((2 * mu_x.0 * mu_y.0) as f32 + c1) * ((2 * sigma_xy.0) as f32 + c2))
            / (((mu_x.0.pow(2) + mu_y.0.pow(2)) as f32 + c1)
                * ((sigma_x.0 + sigma_y.0) as f32 + c2));

        let ssim_g = (((2 * mu_x.1 * mu_y.1) as f32 + c1) * ((2 * sigma_xy.1) as f32 + c2))
            / (((mu_x.1.pow(2) + mu_y.1.pow(2)) as f32 + c1)
                * ((sigma_x.1 + sigma_y.1) as f32 + c2));

        let ssim_b = (((2 * mu_x.2 * mu_y.2) as f32 + c1) * ((2 * sigma_xy.2) as f32 + c2))
            / (((mu_x.2.pow(2) + mu_y.2.pow(2)) as f32 + c1)
                * ((sigma_x.2 + sigma_y.2) as f32 + c2));

        (ssim_r + ssim_g + ssim_b) / 3.0
    }
}

#[derive(Clone, Hash, Copy)]
struct Point {
    pub x: u32,
    pub y: u32,
}

impl Point {
    pub fn new(x: u32, y: u32) -> Self {
        Point { x, y }
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Point) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl Eq for Point {}
