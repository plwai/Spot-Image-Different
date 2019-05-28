use std::hash::{Hash, Hasher};

#[derive(Clone, Copy)]
pub struct Point {
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

    pub fn calculate_distance(&self, point_2: &Point) -> f32 {
        let (x, y) = (
            (self.x as f32 - point_2.x as f32),
            (self.y as f32 - point_2.y as f32),
        );

        x.hypot(y)
    }

    pub fn calculate_difference(&self, point_2: &Point) -> Self {
        let x = self.x - point_2.x;
        let y = self.y - point_2.y;

        Point {
            x_id: x as u32,
            y_id: y as u32,
            x,
            y,
        }
    }

    pub fn cross_product(&self, point_2: &Point) -> f32 {
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
