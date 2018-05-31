extern crate cgmath;
extern crate image;
extern crate rayon;
extern crate chrono;

mod tracer;
use tracer::scene::*;
use tracer::shapes::*;

use std::fs::File;

use cgmath::prelude::*;
use cgmath::Point3;
use cgmath::Vector2;
use cgmath::Vector3;

use rayon::prelude::*;

use chrono::prelude::*;

const GAMMA: f32 = 2.2;

fn gamma_encode(linear: f32) -> f32 {
    linear.powf(1.0 / GAMMA)
}

fn gamma_decode(encoded: f32) -> f32 {
    encoded.powf(GAMMA)
}


fn image_correction(pixels: Vec<f32>) -> Vec<u8> {
    let max_value = pixels.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    pixels.iter().enumerate().map(|(i, v)| {
        if i % 4 != 3 {
            (gamma_encode(v / max_value) * 255.0) as u8
        }
        else {
            255
        }
    }).collect()
}

fn main() {
    let scene = Scene {
        materials: vec![
            Material {
                surface_color: Color::new(1.0, 0.32, 0.36),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                emission_color: Color::zero()
            },
            Material {
                surface_color: Color::new(0.90, 0.76, 0.46),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                emission_color: Color::zero()
            },
            Material {
                surface_color: Color::new(0.32, 0.36, 0.90),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                emission_color: Color::zero()
            },
            Material {
                surface_color: Color::new(1.0, 1.0, 0.3),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(1.0, 1.0, 0.3)
            },
            Material {
                surface_color: Color::new(0.0, 0.0, 0.0),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(1.5, 1.5, 1.5)
            },
            Material {
                surface_color: Color::new(0.8, 0.8, 0.8),
                reflectivity: 0.7,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(0.0, 0.0, 0.0)
            }
        ],
        models: vec![
            Model::cube(2)
        ],
        spheres: vec![
            Sphere {
                origin: Point3::new(-1.0, -2.0, -3.0),
                radius: 2.0,
                mat_id: 0
            },
            Sphere {
                origin: Point3::new(2.0, -2.0, -6.0),
                radius: 2.0,
                mat_id: 1
            },
            Sphere {
                origin: Point3::new(-4.0, -2.0, -8.0),
                radius: 2.0,
                mat_id: 3
            },
            Sphere {
                origin: Point3::new(0.0, 0.0, -10020.0),
                radius: 10000.0,
                mat_id: 4
            },
            Sphere {
                origin: Point3::new(0.0, -10004.0, 0.0),
                radius: 10000.0,
                mat_id: 5
            }
        ],
        camera_pos: Point3::new(1.2, 0.7, 10.0)
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
    let mut pixels = vec![0.0; (image_width * image_height * 4) as usize];

    let time = Local::now();

    if multithreading {
        let bands: Vec<(usize, &mut [f32])> = pixels
            .chunks_mut((image_width * 4) as usize)
            .enumerate()
            .collect();

        bands.into_par_iter().for_each(|(i, band)| {
            let band_upper_left = (0, i as u32);
            let band_lower_right = (image_width, i as u32 + 1);
            scene.render(band, (image_width, image_height), band_upper_left, band_lower_right);
        });
    }
    else {
        let time = Local::now();
        scene.render(&mut pixels, (image_width, image_height), (0, 0), (image_width, image_height));
    }

    println!("Elapsed time: {}ms", Local::now().signed_duration_since(time).num_milliseconds());

    let image_data = image_correction(pixels);

    image::save_buffer("result.png", &image_data, 1920, 1080, image::RGBA(8)).unwrap();
}
