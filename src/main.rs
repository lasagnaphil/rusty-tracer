extern crate cgmath;
extern crate image;
extern crate rayon;
extern crate chrono;
extern crate tobj;
#[macro_use] extern crate itertools;
extern crate pbr;
#[macro_use] extern crate lazy_static;
extern crate clap;

mod tracer;
use tracer::scene::*;
use tracer::shapes::*;
use tracer::loader::*;

use std::fs::File;
use std::sync::{Arc, Mutex};

use cgmath::prelude::*;
use cgmath::Point3;
use cgmath::Vector2;
use cgmath::Vector3;

use rayon::prelude::*;

use chrono::prelude::*;

use pbr::ProgressBar;

const GAMMA: f32 = 2.2;

fn gamma_encode(linear: f32) -> f32 {
    linear.powf(1.0 / GAMMA)
}

fn gamma_decode(encoded: f32) -> f32 {
    encoded.powf(GAMMA)
}

fn image_correction(pixels: Vec<f32>) -> Vec<u8> {
    let max_value = pixels.iter().max_by(|a, b| {
        a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal)
    }).unwrap();
    pixels.iter().enumerate().map(|(i, v)| {
        if i % 4 != 3 {
            (gamma_encode(v / max_value) * 255.0) as u8
        }
        else {
            255
        }
    }).collect()
}

use clap::{Arg, App};

fn main() {
    let mut scene = Scene {
        materials: vec![
            Material {
                surface_color: Color::new(0.9, 0.32, 0.36),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::zero()
            },
            Material {
                surface_color: Color::new(0.2, 0.76, 0.1),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                emission_color: Color::zero()
            },
            Material {
                surface_color: Color::new(0.9, 0.9, 0.2),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                emission_color: Color::new(0.0, 0.0, 0.0)
            },
            Material {
                surface_color: Color::new(0.1, 0.2, 0.9),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(0.0, 0.0, 0.0)
            },
            Material {
                surface_color: Color::new(0.8, 0.8, 0.8),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(0.0, 0.0, 0.0)
            },
            Material {
                surface_color: Color::new(0.0, 0.0, 0.0),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(1.0, 1.0, 1.0)
            }
        ],
        meshes: vec![
        ],
        spheres: vec![
            /*
            Sphere {
                origin: Point3::new(-1.0, -2.0, -3.0),
                radius: 2.0,
                mat_id: 0
            },
            Sphere {
                origin: Point3::new(3.0, -2.0, -3.0),
                radius: 2.0,
                mat_id: 1
            },
            Sphere {
                origin: Point3::new(2.0, -2.0, -6.0),
                radius: 2.0,
                mat_id: 2
            },
            Sphere {
                origin: Point3::new(7.0, -2.0, -8.0),
                radius: 2.0,
                mat_id: 3
            },
            Sphere {
                origin: Point3::new(0.0, -10004.0, 0.0),
                radius: 10000.0,
                mat_id: 4
            },
            */
            /*
            Sphere {
                origin: Point3::new(0.0, 0.0, -10020.0),
                radius: 10000.0,
                mat_id: 5
            },
            Sphere {
                origin: Point3::new(0.0, 0.0, -10020.0),
                radius: 10000.0,
                mat_id: 5
            },
            */
        ],
        point_lights: vec![
            PointLight {
                pos: Point3f::new(10.0, 10.0, 10.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            PointLight {
                pos: Point3f::new(-10.0, 10.0, 10.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            /*
            PointLight {
                pos: Point3f::new(0.0, 0.0, 1000.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            PointLight {
                pos: Point3f::new(0.0, 0.0, -1000.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            PointLight {
                pos: Point3f::new(1000.0, 0.0, 0.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            PointLight {
                pos: Point3f::new(-1000.0, 0.0, 0.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            PointLight {
                pos: Point3f::new(0.0, 1000.0, 0.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            PointLight {
                pos: Point3f::new(0.0, -1000.0, 0.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            */
        ],
        camera_pos: Point3::new(0.0, 2.0, 10.0)
    };

    let matches = App::new("rusty-tracer")
        .version("1.0")
        .author("Philsik Chang <lasagnaphil@snu.ac.kr>")
        .about("A ray tracer written in rust")
        .arg(Arg::with_name("single-thread")
            .short("s").long("single-thread"))
        .arg(Arg::with_name("INPUT").index(1))
        .get_matches();

    match matches.value_of("INPUT") {
        Some(filename) => {
            let mut obj_meshes = obj_to_meshes(filename);
            println!("Model {} successfully loaded.", filename);
            for mesh in obj_meshes {
                scene.add_mesh(mesh);
            }
        }
        None => {}
    }
    let multithreading = !matches.is_present("single-thread");

    let mut cube_mesh = Mesh::cube(2);
    let transform = Matrix4f::from_angle_y(cgmath::Deg(30.0));
    let transform = transform * Matrix4f::from_angle_z(cgmath::Deg(30.0));
    let transform = transform * Matrix4f::from_translation(Vector3::new(0.5, -1.0, 0.0));
    cube_mesh.transform(transform);
    // scene.add_mesh(cube_mesh);

    let image_width: u32 = 1280;
    let image_height: u32 = 720;
    let mut pixels = vec![0.0; (image_width * image_height * 4) as usize];

    let time = Local::now();

    if multithreading {
        println!("Starting ray tracer with multiple threads.");
        let bands: Vec<(usize, &mut [f32])> = pixels
            .chunks_mut((image_width * 4) as usize)
            .enumerate()
            .collect();

        let mut pb = Arc::new(Mutex::new(ProgressBar::new(bands.len() as u64)));
        pb.lock().unwrap().format("╢▌▌░╟");

        bands.into_par_iter().for_each(|(i, band)| {
            let pb = pb.clone();
            pb.lock().unwrap().inc();
            let band_upper_left = (0, i as u32);
            let band_lower_right = (image_width, i as u32 + 1);
            scene.render(band, (image_width, image_height), band_upper_left, band_lower_right);
        });

        println!("Ray tracing complete!")
    }
    else {
        println!("Starting ray tracer with single thread.");

        let time = Local::now();
        scene.render(&mut pixels, (image_width, image_height), (0, 0), (image_width, image_height));

        println!("Ray tracing complete!")
    }

    println!("Elapsed time: {}ms", Local::now().signed_duration_since(time).num_milliseconds());

    println!("Performing postprocessing.");
    let image_data = image_correction(pixels);

    println!("Saving image.");
    image::save_buffer("result.png", &image_data, image_width, image_height, image::RGBA(8)).unwrap();

    println!("Rendered image saved at result.png");
}
