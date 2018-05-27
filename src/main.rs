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
type Color = Vector3<f32>;

const GAMMA: f32 = 2.2;

fn gamma_encode(linear: f32) -> f32 {
    linear.powf(1.0 / GAMMA)
}

fn gamma_decode(encoded: f32) -> f32 {
    encoded.powf(GAMMA)
}

pub fn color_to_rgba(color: &Color) -> Rgba<u8> {
    Rgba::from_channels((gamma_encode(color.x.min(1.0)) * 255.0) as u8,
                        (gamma_encode(color.y.min(1.0)) * 255.0) as u8,
                        (gamma_encode(color.z.min(1.0)) * 255.0) as u8,
                        255)
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
    surface_color: Color,
    reflectivity: f32,
    transparency: f32,
    refractive_index: f32,
    emission_color: Color
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

fn trace(image: &mut DynamicImage, scene: &Scene, ray: &Ray, depth: u32) -> Color {

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
        if (sphere.transparency > 0.0 || sphere.reflectivity > 0.0) && depth < MAX_RAY_DEPTH {
            let n = if incident_angle > 0.0 {
                1.0 / sphere.refractive_index
            } else {
                sphere.refractive_index
            };
            let r0 = ((n - 1.0) / (n + 1.0)).powi(2);
            let fresnel = r0 + (1.0 - r0) * (1.0 - incident_angle.abs()).powi(5);

            let reflect_dir = reflect(ray.dir, hit_normal);
            let reflection_ray = Ray::new(hit_pos + hit_normal * 1e-4, reflect_dir);
            let reflection_color = trace(image, scene, &reflection_ray, depth + 1);

            let refraction_color = if sphere.transparency > 0.0 {
                let refract_dir = refract(ray.dir, hit_normal, n);
                let refraction_ray = Ray::new(hit_pos - hit_normal * 1e-4, refract_dir);
                trace(image, scene, &refraction_ray, depth + 1)
            } else { Color::zero() };

            sphere.emission_color + sphere.surface_color.mul_element_wise(
                reflection_color * fresnel + refraction_color * (1.0 - fresnel) * sphere.transparency
            )
        } else {
            let mut surface_color = Color::zero();
            for light_sphere in &scene.spheres {
                if light_sphere.emission_color != Color::zero() {
                    let shadow_ray = Ray::new(hit_pos, (light_sphere.origin - hit_pos).normalize());
                    let is_shadow = scene.spheres.iter().any(|other_sphere| {
                        intersect_sphere(other_sphere, &shadow_ray).is_some()
                    });
                    if !is_shadow {
                        let shadow_angle = hit_normal.dot(shadow_ray.dir);
                        if shadow_angle > 0.0 {
                            surface_color += shadow_angle *
                                sphere.surface_color.mul_element_wise(light_sphere.emission_color);
                        }
                    };
                }
            }
            sphere.emission_color + surface_color
        }
    }
    else {
        Vector3f::zero()
    }
}

fn render(scene: &Scene, image_width: u32, image_height: u32) -> DynamicImage {
    let mut image = DynamicImage::new_rgb8(image_width, image_height);

    let fov = 90.0f32;
    let tangent = (fov / 2.0f32).to_radians().tan();
    let aspect_ratio = image_height as f32 / image_width as f32;

    for j in 0..image_height {
        for i in 0..image_width {
            let prim_ray_dir = Vector3::new(
                ((2.0 * i as f32 / image_width as f32) - 1.0) * tangent,
                -((2.0 * j as f32 / image_width as f32) - aspect_ratio) * tangent,
                -1.0).normalize();

            let prim_ray = Ray::new(scene.camera_pos, prim_ray_dir);

            let color = trace(&mut image, &scene, &prim_ray, 0);
            image.put_pixel(i, j, color_to_rgba(&color));
        }
    }
    return image;
}

fn main() {
    let scene = Scene {
        spheres: vec![
            Sphere {
                origin: Point3::new(-1.0, -2.0, -3.0),
                radius: 2.0,
                surface_color: Color::new(1.0, 0.32, 0.36),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                emission_color: Color::zero()
            },
            Sphere {
                origin: Point3::new(2.0, -2.0, -6.0),
                radius: 2.0,
                surface_color: Color::new(0.90, 0.76, 0.46),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                emission_color: Color::zero()
            },
            Sphere {
                origin: Point3::new(-4.0, -2.0, -8.0),
                radius: 2.0,
                surface_color: Color::new(1.0, 1.0, 0.3),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(3.0, 3.0, 0.9)
            },
            Sphere {
                origin: Point3::new(0.0, 0.0, -10020.0),
                radius: 10000.0,
                surface_color: Color::new(0.0, 0.0, 0.0),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(0.5, 0.5, 0.5)
            },
            Sphere {
                origin: Point3::new(0.0, -10004.0, 0.0),
                radius: 10000.0,
                surface_color: Color::new(0.8, 0.8, 0.8),
                reflectivity: 0.7,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(0.0, 0.0, 0.0)
            }
        ],
        camera_pos: Point3::new(0.0, 0.0, 10.0)
    };
    let image = render(&scene, 1920, 1080);
    image.save("result.png").unwrap();
}
