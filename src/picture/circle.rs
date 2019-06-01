use rand::{seq::SliceRandom, thread_rng};
use std::collections::HashMap;

use super::Point;

#[derive(Clone, Copy)]
pub struct Circle {
    pub mid_point: Point,
    pub radius: f32,
}

impl Circle {
    pub fn new(mid_point: Point, radius: f32) -> Self {
        Circle { mid_point, radius }
    }

    pub fn contains_point(&self, point: &Point) -> bool {
        self.mid_point.distance(point) <= self.radius * (1.0 + 1e-14)
    }
}

// Smallest Circle Problem
pub fn compute_smallest_circle(points_group: HashMap<i32, Vec<Point>>) -> Vec<Option<Circle>> {
    let circle = points_group
        .iter()
        .map(|(_group, points)| {
            let shuffled_points = &mut points.clone()[..];
            let mut rng = thread_rng();
            shuffled_points.shuffle(&mut rng);

            let shuffled_points = shuffled_points.to_vec();

            let mut circle: Option<Circle> = None;

            for (index, point) in shuffled_points.iter().enumerate() {
                if circle.is_some() && circle.unwrap().contains_point(point) {
                    continue;
                }

                circle = construct_one_point_circle_enclosing(&shuffled_points, point, index + 1);
            }

            circle
        })
        .collect();

    circle
}

fn construct_diameter(point_1: &Point, point_2: &Point) -> Circle {
    let x = (point_1.x + point_2.x) / 2.0;
    let y = (point_1.y + point_2.y) / 2.0;

    let new_point = Point::new(x, y);

    let radius = new_point
        .distance(point_1)
        .max(new_point.distance(point_2));

    Circle::new(new_point, radius)
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
                circle = construct_diameter(point_1, point_2);
            } else {
                circle = construct_two_point_circle_enclosing(
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
    let circle = construct_diameter(point_1, point_2);
    let mut left: Option<Circle> = None;
    let mut right: Option<Circle> = None;

    let point_3 = point_1.difference(point_2);

    for index in 0..end {
        let point_r = &shuffled_points[index];

        if circle.contains_point(point_r) {
            continue;
        }

        let cross_result = point_3.cross_product(&point_r.difference(point_1));
        let circle_2 = construct_circumcircle(point_1, point_2, point_r);

        match circle_2 {
            Some(c) => {
                if cross_result > 0.0
                    && (left.is_none()
                        || point_3.cross_product(&c.mid_point.difference(point_1))
                            > point_3.cross_product(
                                &left.unwrap().mid_point.difference(point_1),
                            ))
                {
                    left = Some(c);
                } else if cross_result < 0.0
                    && (right.is_none()
                        || point_3.cross_product(&c.mid_point.difference(point_1))
                            > point_3.cross_product(
                                &right.unwrap().mid_point.difference(point_1),
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
        .distance(a)
        .max(p.distance(b))
        .max(p.distance(c));

    Some(Circle {
        mid_point: p,
        radius: r,
    })
}
