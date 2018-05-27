extern crate cgmath;
extern crate image;
extern crate rayon;
extern crate chrono;

use std::f32;
use std::fs::File;

use cgmath::prelude::*;
use cgmath::Point3;
use cgmath::Vector3;

use image::ColorType;
use image::png::PNGEncoder;

use rayon::prelude::*;

use chrono::prelude::*;

type Point3f = Point3<f32>;
type Vector3f = Vector3<f32>;
type Color = Vector3<f32>;

const GAMMA: f32 = 2.2;

fn gamma_encode(linear: f32) -> f32 {
    linear.powf(1.0 / GAMMA)
}

fn gamma_decode(encoded: f32) -> f32 {
    encoded.powf(GAMMA)
}

pub fn color_to_rgba(color: &Color) -> [u8; 4] {
    [
        (gamma_encode(color.x.min(1.0)) * 255.0) as u8,
        (gamma_encode(color.y.min(1.0)) * 255.0) as u8,
        (gamma_encode(color.z.min(1.0)) * 255.0) as u8,
        255
    ]
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
pub struct Material {
    surface_color: Color,
    reflectivity: f32,
    transparency: f32,
    refractive_index: f32,
    emission_color: Color
}

#[derive(Clone)]
pub struct Sphere {
    origin: Point3f,
    radius: f32,
    mat: Material
}

pub struct Scene {
    camera_pos: Point3f,
    spheres: Vec<Sphere>
}

fn intersect_sphere(sphere: &Sphere, ray: &Ray) -> Option<(f32, f32)> {
    let disp = sphere.origin - ray.origin;
    let ip = ray.dir.dot(disp);
    if ip < 0.0 { return None; }
    let discriminant = ip * ip - disp.magnitude2() + sphere.radius * sphere.radius;
    if discriminant >= 0.0 {
        Some((ip - discriminant.sqrt(), ip + discriminant.sqrt()))
    }
    else {
        None
    }
}

fn reflect(dir: Vector3f, normal: Vector3f) -> Vector3f {
    (dir - 2.0 * (dir.dot(normal)) * normal).normalize()
}

fn refract(dir: Vector3f, normal: Vector3f, n: f32) -> Vector3f {
    let cos_theta_i = -dir.dot(normal);
    let cos_theta_r = (1.0 - n*n*(1.0 - cos_theta_i*cos_theta_i)).sqrt() as f32;
    ((n * cos_theta_i - cos_theta_r) * normal + n * dir).normalize()
}

const MAX_RAY_DEPTH: u32 = 5;

fn trace(scene: &Scene, ray: &Ray, depth: u32) -> Color {

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

    if let Some(sphere) = closest_sphere {
        let hit_pos = ray.origin + tnear * ray.dir;
        let hit_normal = (hit_pos - sphere.origin).normalize();
        let incident_angle = -ray.dir.dot(hit_normal);
        let hit_normal = if incident_angle > 0.0 { hit_normal } else { -hit_normal };
        if (sphere.mat.transparency > 0.0 || sphere.mat.reflectivity > 0.0) && depth < MAX_RAY_DEPTH {
            let n = if incident_angle > 0.0 {
                1.0 / sphere.mat.refractive_index
            } else {
                sphere.mat.refractive_index
            };
            let r0 = ((n - 1.0) / (n + 1.0)).powi(2);
            let fresnel = r0 + (1.0 - r0) * (1.0 - incident_angle.abs()).powi(5);

            let reflect_dir = reflect(ray.dir, hit_normal);
            let reflection_ray = Ray::new(hit_pos + hit_normal * 1e-4, reflect_dir);
            let reflection_color = trace(scene, &reflection_ray, depth + 1);

            let refraction_color = if sphere.mat.transparency > 0.0 {
                let refract_dir = refract(ray.dir, hit_normal, n);
                let refraction_ray = Ray::new(hit_pos - hit_normal * 1e-4, refract_dir);
                trace(scene, &refraction_ray, depth + 1)
            } else { Color::zero() };

            sphere.mat.emission_color + sphere.mat.surface_color.mul_element_wise(
                reflection_color * fresnel + refraction_color * (1.0 - fresnel) * sphere.mat.transparency
            )
        } else {
            let mut surface_color = Color::zero();
            for light_sphere in &scene.spheres {
                if light_sphere.mat.emission_color != Color::zero() {
                    let shadow_ray = Ray::new(hit_pos, (light_sphere.origin - hit_pos).normalize());
                    let is_shadow = scene.spheres.iter().any(|other_sphere| {
                        intersect_sphere(other_sphere, &shadow_ray).is_some()
                    });
                    if !is_shadow {
                        let shadow_angle = hit_normal.dot(shadow_ray.dir);
                        if shadow_angle > 0.0 {
                            surface_color += shadow_angle *
                                sphere.mat.surface_color.mul_element_wise(light_sphere.mat.emission_color);
                        }
                    };
                }
            }
            sphere.mat.emission_color + surface_color
        }
    }
    else {
        Vector3f::zero()
    }
}

fn render(pixels: &mut [u8], scene: &Scene, bounds: (u32, u32), upper_left: (u32, u32), lower_right: (u32, u32)) {
    let fov = 60.0f32;
    let tangent = (fov / 2.0f32).to_radians().tan();
    let aspect_ratio = bounds.1 as f32 / bounds.0 as f32;

    let iter_width = lower_right.0 - upper_left.0;
    let iter_height = lower_right.1 - upper_left.1;
    for j in 0..iter_height {
        for i in 0..iter_width {
            let prim_ray_dir = Vector3::new(
                ((2.0 * (upper_left.0 + i) as f32 / bounds.0 as f32) - 1.0) * tangent,
                -((2.0 * (upper_left.1 + j) as f32 / bounds.0 as f32) - aspect_ratio) * tangent,
                -1.0).normalize();

            let prim_ray = Ray::new(scene.camera_pos, prim_ray_dir);

            let color = trace(&scene, &prim_ray, 0);
            let color = color_to_rgba(&color);
            pixels[(4*(j * iter_width + i)) as usize] = color[0];
            pixels[(4*(j * iter_width + i) + 1) as usize] = color[1];
            pixels[(4*(j * iter_width + i) + 2) as usize] = color[2];
            pixels[(4*(j * iter_width + i) + 3) as usize] = color[3];
        }
    }
}

fn main() {
    let scene = Scene {
        spheres: vec![
            Sphere {
                origin: Point3::new(-1.0, -2.0, -3.0),
                radius: 2.0,
                mat: Material {
                    surface_color: Color::new(1.0, 0.32, 0.36),
                    reflectivity: 1.0,
                    transparency: 0.5,
                    refractive_index: 1.1,
                    emission_color: Color::zero()
                }
            },
            Sphere {
                origin: Point3::new(2.0, -2.0, -6.0),
                radius: 2.0,
                mat: Material {
                    surface_color: Color::new(0.90, 0.76, 0.46),
                    reflectivity: 1.0,
                    transparency: 0.5,
                    refractive_index: 1.1,
                    emission_color: Color::zero()
                }
            },
            Sphere {
                origin: Point3::new(-4.0, -2.0, -8.0),
                radius: 2.0,
                mat: Material {
                    surface_color: Color::new(1.0, 1.0, 0.3),
                    reflectivity: 0.0,
                    transparency: 0.0,
                    refractive_index: 1.1,
                    emission_color: Color::new(3.0, 3.0, 0.9)
                }
            },
            Sphere {
                origin: Point3::new(0.0, 0.0, -10020.0),
                radius: 10000.0,
                mat: Material {
                    surface_color: Color::new(0.0, 0.0, 0.0),
                    reflectivity: 0.0,
                    transparency: 0.0,
                    refractive_index: 1.1,
                    emission_color: Color::new(0.5, 0.5, 0.5)
                }
            },
            Sphere {
                origin: Point3::new(0.0, -10004.0, 0.0),
                radius: 10000.0,
                mat: Material {
                    surface_color: Color::new(0.8, 0.8, 0.8),
                    reflectivity: 0.7,
                    transparency: 0.0,
                    refractive_index: 1.1,
                    emission_color: Color::new(0.0, 0.0, 0.0)
                }
            }
        ],
        camera_pos: Point3::new(0.0, 0.0, 10.0)
    };

    let args: Vec<String> = std::env::args().collect();
    let multithreading = match args.len() {
        1 => {
            true
        }
        2 => {
            if args[1] == "--single-thread" {
                false
            } else {
                eprintln!("Invalid argument!");
                std::process::exit(1);
            }
        }
        _ => {
            eprintln!("Invalid argument!");
            std::process::exit(1);
        }
    };

    let image_width: u32 = 1920;
    let image_height: u32 = 1080;
    let mut pixels = vec![0; (image_width * image_height * 4) as usize];

    if multithreading {
        let time = Local::now();
        let bands: Vec<(usize, &mut [u8])> = pixels
            .chunks_mut((image_width * 4) as usize)
            .enumerate()
            .collect();

        bands.into_par_iter().for_each(|(i, band)| {
            let band_upper_left = (0, i as u32);
            let band_lower_right = (image_width, i as u32 + 1);
            render(band, &scene, (image_width, image_height), band_upper_left, band_lower_right);
        });

        println!("Elapsed time: {}ms", Local::now().signed_duration_since(time).num_milliseconds());
    }
    else {
        let time = Local::now();
        render(&mut pixels, &scene, (image_width, image_height), (0, 0), (image_width, image_height));

        println!("Elapsed time: {}ms", Local::now().signed_duration_since(time).num_milliseconds());
    }

    image::save_buffer("result.png", &pixels, 1920, 1080, image::RGBA(8)).unwrap();
}
