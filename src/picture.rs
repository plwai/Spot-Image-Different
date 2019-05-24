use image;
use image::*;

use rand::{seq::SliceRandom, thread_rng};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

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
        let points =
            Picture::get_diff_points(&mut self.img, &mut rhs_picture.img, (10, 10), 0.4, 0.9);

        println!("");

        let point_groups = Picture::dbscan(&points, 20.0, 1);
        let result = point_groups
            .iter()
            .fold(HashMap::new(), |mut acc, (point, group)| {
                let group = group.unwrap();
                let points = acc.entry(group).or_insert(vec![]);
                points.push(*point);
                acc
            });

        let mut count = 0;
        let mut group_num = 1;

        for (group, points) in result.iter() {
            println!("");
            println!("");
            println!("Group Num {}: Cluster Group: {}", group_num, group);
            for p in points.iter() {
                print!("{} {}, ", p.x, p.y);
                count = count + 1;
            }
            group_num = group_num + 1;
        }
        println!("");
        println!("");
        println!("{}", count);

        let circles = Self::compute_smallest_circle(result);

        for circling in circles {
            match circling {
                Some(c) => {
                    println!(
                        "mid point: {} {}, radius: {}",
                        c.mid_point.x, c.mid_point.y, c.radius
                    );
                    rhs_picture.draw_circle(
                        c.mid_point.x as i32,
                        c.mid_point.y as i32,
                        c.radius as i32 + 13,
                        c.radius as i32 + 17,
                        (255, 0, 0),
                    );
                }
                None => println!("NO CIRCLING! PPP"),
            }
        }

        rhs_picture.img.save("Answer.jpg").unwrap();
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
                    //println!("{} {} {}", current_origin.0, current_origin.1, ssim);

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
                                current_origin.0 as f32 + point.x,
                                current_origin.1 as f32 + point.y,
                            ));
                        }
                    } else {
                        local_points
                            .push(Point::new(current_origin.0 as f32, current_origin.1 as f32));
                    }
                }
            }
        }

        local_points
    }

    fn dbscan(
        points: &Vec<Point>,
        eps: f32,
        min_pts: usize,
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
                .filter(|p| p.calculate_distance(&point) <= eps)
                .collect::<Vec<_>>();

            if neighbors.len() < min_pts {
                // -1 as noise
                points_map.insert(*point, Some(-1));
                continue;
            }

            cluster_id = cluster_id + 1;
            points_map.insert(*point, Some(cluster_id));

            let mut start = 0;

            loop {
                for index in start..neighbors.len() {
                    start = index;
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
                        .filter(|p| p.calculate_distance(sp) <= eps && !neighbors.contains(p))
                        .collect::<Vec<_>>();

                    if seed_neighbors.len() >= min_pts {
                        neighbors.extend(seed_neighbors);
                    }
                }

                if start == neighbors.len() - 1 {
                    break;
                }
            }
        }

        points_map
    }

    fn construct_diameter(point_1: &Point, point_2: &Point) -> Circle {
        let x = (point_1.x + point_2.x) / 2.0;
        let y = (point_1.y + point_2.y) / 2.0;

        let new_point = Point::new(x, y);

        let radius = new_point
            .calculate_distance(point_1)
            .max(new_point.calculate_distance(point_2));

        Circle::new(new_point, radius)
    }

    // Smallest Circle Problem
    fn compute_smallest_circle(points_group: HashMap<i32, Vec<Point>>) -> Vec<Option<Circle>> {
        let circle = points_group
            .iter()
            .map(|(group, points)| {
                let shuffled_points = &mut points.clone()[..];
                let mut rng = thread_rng();
                shuffled_points.shuffle(&mut rng);

                let shuffled_points = shuffled_points.to_vec();

                let mut circle: Option<Circle> = None;

                for (index, point) in shuffled_points.iter().enumerate() {
                    if circle.is_some() && circle.unwrap().contains_point(point) {
                        continue;
                    }

                    circle = Self::construct_one_point_circle_enclosing(
                        &shuffled_points,
                        point,
                        index + 1,
                    );
                }

                circle
            })
            .collect();

        circle
    }

    fn construct_one_point_circle_enclosing(
        shuffled_points: &Vec<Point>,
        point_1: &Point,
        end: usize,
    ) -> Option<Circle> {
        let mut circle = Circle::new(*point_1, 0.0);

        for index in 0..end {
            let point_2 = &shuffled_points[index];

            if !circle.contains_point(point_2) {
                if circle.radius == 0.0 {
                    circle = Self::construct_diameter(point_1, point_2);
                } else {
                    circle = Self::construct_two_point_circle_enclosing(
                        shuffled_points,
                        point_1,
                        point_2,
                        index + 1,
                    )
                    .unwrap();
                }
            }
        }

        Some(circle)
    }

    fn construct_two_point_circle_enclosing(
        shuffled_points: &Vec<Point>,
        point_1: &Point,
        point_2: &Point,
        end: usize,
    ) -> Option<Circle> {
        let circle = Self::construct_diameter(point_1, point_2);
        let mut left: Option<Circle> = None;
        let mut right: Option<Circle> = None;

        let point_3 = point_1.calculate_difference(point_2);

        for index in 0..end {
            let point_r = &shuffled_points[index];

            if circle.contains_point(point_r) {
                continue;
            }

            let cross_result = point_3.cross_product(&point_r.calculate_difference(point_1));
            let circle_2 = Self::construct_circumcircle(point_1, point_2, point_r);

            match circle_2 {
                Some(c) => {
                    if cross_result > 0.0
                        && (left.is_none()
                            || point_3.cross_product(&c.mid_point.calculate_difference(point_1))
                                > point_3.cross_product(
                                    &left.unwrap().mid_point.calculate_difference(point_1),
                                ))
                    {
                        left = Some(c);
                    } else if cross_result < 0.0
                        && (right.is_none()
                            || point_3.cross_product(&c.mid_point.calculate_difference(point_1))
                                > point_3.cross_product(
                                    &right.unwrap().mid_point.calculate_difference(point_1),
                                ))
                    {
                        right = Some(c);
                    }
                }
                None => continue,
            }
        }

        if left.is_none() && right.is_none() {
            return Some(circle);
        } else if left.is_none() {
            return right;
        } else if right.is_none() {
            return left;
        } else {
            if left.unwrap().radius <= right.unwrap().radius {
                return left;
            } else {
                return right;
            };
        }
    }

    fn construct_circumcircle(a: &Point, b: &Point, c: &Point) -> Option<Circle> {
        let ox = (a.x.min(b.x).min(c.x) + a.x.min(b.x).max(c.x)) / 2.0;
        let oy = (a.y.min(b.y).min(c.y) + a.y.min(b.y).max(c.y)) / 2.0;
        let (ax, ay) = (a.x - ox, a.y - oy);
        let (bx, by) = (b.x - ox, b.y - oy);
        let (cx, cy) = (c.x - ox, c.y - oy);
        let d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0;

        if d == 0.0 {
            return None;
        }

        let x = ((ax * ax + ay * ay) * (by - cy)
            + (bx * bx + by * by) * (cy - ay)
            + (cx * cx + cy * cy) * (ay - by))
            / d;
        let y = ((ax * ax + ay * ay) * (cx - bx)
            + (bx * bx + by * by) * (ax - cx)
            + (cx * cx + cy * cy) * (bx - ax))
            / d;;

        let p = Point::new(ox + x, oy + y);
        let r = p
            .calculate_distance(a)
            .max(p.calculate_distance(b))
            .max(p.calculate_distance(c));

        Some(Circle {
            mid_point: p,
            radius: r,
        })
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

#[derive(Clone, Copy)]
struct Point {
    pub x_id: u32,
    pub y_id: u32,
    pub x: f32,
    pub y: f32,
}

impl Point {
    pub fn new(x: f32, y: f32) -> Self {
        Point {
            x_id: x as u32,
            y_id: y as u32,
            x,
            y,
        }
    }

    fn calculate_distance(&self, point_2: &Point) -> f32 {
        let (x, y) = (
            (self.x as f32 - point_2.x as f32),
            (self.y as f32 - point_2.y as f32),
        );

        x.hypot(y)
    }

    fn calculate_difference(&self, point_2: &Point) -> Self {
        let x = self.x - point_2.x;
        let y = self.y - point_2.y;

        Point {
            x_id: x as u32,
            y_id: y as u32,
            x,
            y,
        }
    }

    fn cross_product(&self, point_2: &Point) -> f32 {
        self.x as f32 * point_2.y as f32 - self.y as f32 * point_2.x as f32
    }
}

impl Hash for Point {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.x_id.hash(state);
        self.y_id.hash(state);
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Point) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl Eq for Point {}

#[derive(Clone, Copy)]
struct Circle {
    pub mid_point: Point,
    pub radius: f32,
}

impl Circle {
    pub fn new(mid_point: Point, radius: f32) -> Self {
        Circle { mid_point, radius }
    }

    fn contains_point(&self, point: &Point) -> bool {
        self.mid_point.calculate_distance(point) <= self.radius * (1.0 + 1e-14)
    }
}
