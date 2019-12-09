use std::hash::{Hash, Hasher};
use std::{fmt, fmt::Debug};

#[derive(Clone, Copy)]
pub struct Point {
    pub x_id: i32,
    pub y_id: i32,
    pub x: f32,
    pub y: f32,
}

impl Point {
    pub fn new(x: f32, y: f32) -> Self {
        Point {
            x_id: x as i32,
            y_id: y as i32,
            x,
            y,
        }
    }

    pub fn distance(&self, point_2: &Point) -> f32 {
        let (x, y) = ((self.x - point_2.x), (self.y - point_2.y));

        x.hypot(y)
    }

    pub fn difference(&self, point_2: &Point) -> Self {
        let x = self.x - point_2.x;
        let y = self.y - point_2.y;

        Point {
            x_id: x as i32,
            y_id: y as i32,
            x,
            y,
        }
    }

    pub fn cross_product(&self, point_2: &Point) -> f32 {
        self.x * point_2.y - self.y * point_2.x
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

impl Debug for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Point {{ x: {}, y: {} }}", self.x, self.y)
    }
}

#[cfg(test)]
mod point_test {
    use super::*;

    #[test]
    fn point_creation_test() {
        let point = Point::new(1.0, 2.0);

        assert_eq!(point.x_id, 1);
        assert_eq!(point.y_id, 2);
        assert_eq!(point.x, 1.0);
        assert_eq!(point.y, 2.0);

        let point_expected = Point {
            x_id: 1,
            y_id: 1,
            x: 1.0,
            y: 2.0,
        };

        // Test Eq
        assert!(point == point_expected);
    }

    #[test]
    fn point_distance_normal_test() {
        let point_1 = Point::new(1.0, 2.0);
        let point_2 = Point::new(2.0, 3.0);

        assert_eq!(point_1.distance(&point_2), 1.4142135);
    }

    #[test]
    fn point_distance_same_point_test() {
        let point_1 = Point::new(1.0, 2.0);
        let point_2 = Point::new(1.0, 2.0);

        assert_eq!(point_1.distance(&point_2), 0.0);
    }

    #[test]
    fn point_difference_normal_test() {
        let point_1 = Point::new(1.0, 2.0);
        let point_2 = Point::new(2.0, 3.0);
        let difference_point = Point::new(-1.0, -1.0);

        assert!(point_1.difference(&point_2) == difference_point);
    }

    #[test]
    fn point_difference_same_point_test() {
        let point_1 = Point::new(1.0, 2.0);
        let point_2 = Point::new(1.0, 2.0);
        let difference_point = Point::new(0.0, 0.0);

        assert!(point_1.difference(&point_2) == difference_point);
    }

    #[test]
    fn point_cross_product_normal_test() {
        let point_1 = Point::new(1.0, 2.0);
        let point_2 = Point::new(2.0, 3.0);

        assert_eq!(point_1.cross_product(&point_2), -1.0);
    }

    #[test]
    fn point_cross_product_same_point_test() {
        let point_1 = Point::new(1.0, 2.0);
        let point_2 = Point::new(1.0, 2.0);

        assert_eq!(point_1.cross_product(&point_2), 0.0);
    }
}
