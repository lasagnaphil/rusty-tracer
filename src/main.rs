extern crate cgmath;
extern crate image;

use std::f32;

use cgmath::prelude::*;
use cgmath::Point3;
use cgmath::Vector3;

use image::DynamicImage;
use image::GenericImage;
use image::Rgba;
use image::Pixel;

type Point3f = Point3<f32>;
type Vector3f = Vector3<f32>;

#[derive(Clone)]
pub struct Color {
    r: f32,
    g: f32,
    b: f32
}

const GAMMA: f32 = 2.2;

fn gamma_encode(linear: f32) -> f32 {
    linear.powf(1.0 / GAMMA)
}

fn gamma_decode(encoded: f32) -> f32 {
    encoded.powf(GAMMA)
}

impl Color {
    pub fn new(r: f32, g: f32, b: f32) -> Self {
        Color { r: r, g: g, b: b }
    }
    pub fn to_rgba(&self) -> Rgba<u8> {
        Rgba::from_channels((gamma_encode(self.r) * 255.0) as u8,
                            (gamma_encode(self.g) * 255.0) as u8,
                            (gamma_encode(self.b) * 255.0) as u8,
                            255)
    }
}

pub struct Ray {
    origin: Point3f,
    dir: Vector3f
}

impl Ray {
    pub fn new(origin: Point3f, dir: Vector3f) -> Self {
        Ray { origin: origin, dir: dir }
    }
}

#[derive(Clone)]
pub struct Sphere {
    origin: Point3f,
    radius: f32,
    color: Color,
    is_glass: bool
}

impl Sphere {
    pub fn new(origin: Point3f, radius: f32, color: Color, is_glass: bool) -> Self {
        Sphere { origin, radius, color, is_glass }
    }
}

pub struct Scene {
    camera_pos: Point3f,
    light_pos: Point3f,
    spheres: Vec<Sphere>
}

fn intersect_sphere(sphere: &Sphere, ray: &Ray) -> Option<(f32, f32)> {
    let disp = sphere.origin - ray.origin;
    let ip = ray.dir.dot(disp);
    let discriminant = ip * ip - disp.magnitude2() + sphere.radius * sphere.radius;
    if discriminant >= 0.0 {
        Some((ip - discriminant.sqrt(), ip + discriminant.sqrt()))
    }
    else {
        None
    }
}

fn reflect(dir: Vector3f, normal: Vector3f) -> Vector3f {
    dir - 2.0 * (dir.dot(normal)) * normal
}

fn refract(dir: Vector3f, normal: Vector3f, ni: f32, nr: f32) -> Vector3f {
    let n = ni / nr;
    let cos_theta_i = -dir.dot(normal);
    let cos_theta_r = (1.0 - n*n*(1.0 - cos_theta_i*cos_theta_i)).sqrt();
    (n * cos_theta_i - cos_theta_r) * normal + n * dir
}

const MAX_RAY_DEPTH: u32 = 3;

fn trace(image: &mut DynamicImage, scene: &Scene, ray: &Ray, depth: u32) -> image::Rgba<u8> {

    let mut closest_sphere: Option<Sphere> = None;
    let mut tnear = f32::INFINITY;

    for sphere in &scene.spheres {
        if let Some((tnear_hit, tfar_hit)) = intersect_sphere(sphere, &ray) {
            let t = if tnear_hit < 0.0 { tfar_hit } else { tnear_hit };
            if t < tnear {
                tnear = t;
                closest_sphere = Some(sphere.clone());
            }
        }
    }

    let min_hit_pos = ray.origin + tnear * ray.dir;
    if let Some(sphere) = closest_sphere {
        /*
        if sphere.is_glass && depth < MAX_RAY_DEPTH {
            let reflection_ray = reflect(ray.dir, hit_ray.dir);
        }
        */
        let shadow_ray = Ray::new(min_hit_pos, scene.light_pos - min_hit_pos);
        let is_shadow = scene.spheres.iter().any(|other_sphere| {
            intersect_sphere(other_sphere, &shadow_ray).is_some()
        });
        let is_shadow = false;
        if is_shadow {
            Rgba::from_channels(0, 0, 0, 0)
        }
        else {
            sphere.color.to_rgba()
        }
    }
    else {
        Rgba::from_channels(0, 0, 0, 0)
    }
}

fn render(scene: &Scene, image_width: u32, image_height: u32) -> DynamicImage {
    let depth = 3;

    let camera_pos = Point3::new(0.0, 0.0, 0.0);
    let light_pos = Point3::new(0.0, 2.0, -2.0);

    let mut image = DynamicImage::new_rgb8(image_width, image_height);
    let black = Rgba::from_channels(0, 0, 0, 0);

    let fov = 90.0f32;
    let width_dir_ratio = (fov / 2.0f32).tan();
    let height_dir_ratio = (image_height as f32 / image_width as f32) * width_dir_ratio;

    for j in 0..image_height {
        for i in 0..image_width {
            let prim_ray_dir = Vector3::new(
                i as f32 * width_dir_ratio / (image_width as f32) - width_dir_ratio / 2.0,
                -(j as f32 * height_dir_ratio / (image_height as f32)) + height_dir_ratio / 2.0,
                -1.0).normalize();

            let prim_ray = Ray::new(camera_pos, prim_ray_dir);

            let color = trace(&mut image, &scene, &prim_ray, 0);
            image.put_pixel(i, j, color);
        }
    }
    return image;
}

fn main() {
    let scene = Scene {
        spheres: vec![
            Sphere::new(Point3::new(0.0, 0.0, -5.0),
                        1.0,
                        Color::new(1.0, 0.0, 0.0),
                        false),
            Sphere::new(Point3::new(1.0, 1.0, -3.0),
                        1.0,
                        Color::new(0.0, 1.0, 0.0),
                        true)
        ],
        light_pos: Point3::new(0.0, 0.0, -10.0),
        camera_pos: Point3::new(0.0, 0.0, 0.0)
    };
    let image = render(&scene, 640, 480);
    image.save("result.png").unwrap();
}
