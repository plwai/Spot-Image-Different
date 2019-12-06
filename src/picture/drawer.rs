use crate::picture::Picture;
use image::*;

// https://stackoverflow.com/questions/27755514/circle-with-thickness-drawing-algorithm
pub fn draw_circle(
    picture: &mut Picture,
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
        x_line(picture, xc + xi, xc + xo, yc + y, colour);
        y_line(picture, xc + y, yc + xi, yc + xo, colour);
        x_line(picture, xc - xo, xc - xi, yc + y, colour);
        y_line(picture, xc - y, yc + xi, yc + xo, colour);
        x_line(picture, xc - xo, xc - xi, yc - y, colour);
        y_line(picture, xc - y, yc - xo, yc - xi, colour);
        x_line(picture, xc + xi, xc + xo, yc - y, colour);
        y_line(picture, xc + y, yc - xo, yc - xi, colour);

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

fn draw_pixel(picture: &mut Picture, x: i32, y: i32, colour: (u8, u8, u8)) {
    let mut x = x;
    let mut y = y;

    if x >= picture.img.dimensions().0 as i32 {
        x = picture.img.dimensions().0 as i32 - 1;
    }

    if y >= picture.img.dimensions().1 as i32 {
        y = picture.img.dimensions().1 as i32 - 1;
    }

    if x < 0 {
        x = 0;
    }

    if y < 0 {
        y = 0;
    }

    picture.img.put_pixel(
        x as u32,
        y as u32,
        image::Rgba([colour.0, colour.1, colour.2, 1]),
    );
}

fn x_line(picture: &mut Picture, x1: i32, x2: i32, y: i32, colour: (u8, u8, u8)) {
    let mut x1 = x1;

    while x1 <= x2 {
        draw_pixel(picture, x1, y, colour);
        x1 = x1 + 1;
    }
}

fn y_line(picture: &mut Picture, x: i32, y1: i32, y2: i32, colour: (u8, u8, u8)) {
    let mut y1 = y1;

    while y1 <= y2 {
        draw_pixel(picture, x, y1, colour);
        y1 = y1 + 1;
    }
}
