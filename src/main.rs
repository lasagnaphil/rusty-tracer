extern crate cgmath;
extern crate image;
extern crate rayon;
extern crate chrono;
extern crate tobj;
#[macro_use] extern crate itertools;
extern crate pbr;
#[macro_use] extern crate lazy_static;
extern crate rand;

extern crate clap;

mod tracer;
use tracer::scene::*;
use tracer::shapes::*;
use tracer::loader::*;
use tracer::bvh::*;

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
use rand::prelude::*;

fn main() {
    let mut scene = Scene {
        materials: vec![
            Material {
                surface_color: Color::new(0.9, 0.32, 0.36),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                shininess: 64.0,
                specular_color: Color::new(1.0, 1.0, 1.0),
                ambient_color: Color::zero()
            },
            Material {
                surface_color: Color::new(0.2, 0.76, 0.1),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                shininess: 0.0,
                specular_color: Color::zero(),
                ambient_color: Color::zero()
            },
            Material {
                surface_color: Color::new(0.9, 0.9, 0.2),
                reflectivity: 1.0,
                transparency: 0.5,
                refractive_index: 1.1,
                shininess: 0.0,
                specular_color: Color::zero(),
                ambient_color: Color::zero()
            },
            Material {
                surface_color: Color::new(0.1, 0.2, 0.9),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                shininess: 32.0,
                specular_color: Color::new(0.5, 0.5, 0.5),
                ambient_color: Color::zero()
            },
            Material {
                surface_color: Color::new(0.8, 0.8, 0.8),
                reflectivity: 0.0,
                transparency: 0.0,
                refractive_index: 1.1,
                shininess: 0.0,
                specular_color: Color::zero(),
                ambient_color: Color::zero()
            },
        ],
        meshes: vec![
        ],
        spheres: vec![
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
            /*
            Sphere {
                origin: Point3::new(0.0, 0.0, -10020.0),
                radius: 10000.0,
                mat_id: 5
            },
            Sphere {
                origin: Point3::new(0.0, 0.0, 10020.0),
                radius: 10000.0,
                mat_id: 5
            },
            */
        ],
        point_lights: vec![
            /*
            */
            PointLight {
                pos: Point3f::new(-2.0, 0.0, 1.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            /*
            PointLight {
                pos: Point3f::new(-10.0, 10.0, 10.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
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
        camera_pos: Point3::new(0.0, 2.0, 15.0),
        bvh: None
    };

    let matches = App::new("rusty-tracer")
        .version("1.0")
        .author("Philsik Chang <lasagnaphil@snu.ac.kr>")
        .about("A ray tracer written in rust")
        .arg(Arg::with_name("single-thread")
            .short("s").long("single-thread"))
        .arg(Arg::with_name("bvh")
            .short("b").long("bvh"))
        .arg(Arg::with_name("WIDTH").index(1))
        .arg(Arg::with_name("HEIGHT").index(2))
        .get_matches();

    /*
    match matches.value_of("INPUT") {
        Some(filename) => {
            let mut obj_meshes = obj_to_meshes(filename);
            println!("Model {} successfully loaded.", filename);
            for mut mesh in obj_meshes {
                let transform = Matrix4f::from_angle_x(cgmath::Deg(20.0));
                mesh.transform(transform);
                scene.add_mesh(mesh);
            }
        }
        None => {}
    }
    */

    let multithreading = !matches.is_present("single-thread");
    let use_bvh = matches.is_present("bvh");

    /*
    let mut rng: rand::XorShiftRng = rand::SeedableRng::from_seed(
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    for i in 0..100 {
        let x = rng.gen_range(-10.0, 10.0);
        let y = rng.gen_range(-10.0, 10.0);
        let z = rng.gen_range(-10.0, 10.0);
        let rot_x = rng.gen_range(0.0, 180.0);
        let rot_y = rng.gen_range(0.0, 180.0);
        let rot_z = rng.gen_range(0.0, 180.0);
        let mut cube_mesh = Mesh::cube(0);
        let transform = Matrix4f::from_angle_y(cgmath::Deg(rot_x));
        let transform = transform * Matrix4f::from_angle_y(cgmath::Deg(rot_y));
        let transform = transform * Matrix4f::from_angle_z(cgmath::Deg(rot_z));
        let transform = transform * Matrix4f::from_translation(Vector3::new(x, y, z));
        cube_mesh.transform(transform);
        scene.add_mesh(cube_mesh);
    }
    */

    if use_bvh {
        println!("Building BVH.");
        scene.build_bvh();
        println!("BVH building complete.");
    }

    let image_width: u32 = matches.value_of("WIDTH")
        .unwrap().parse::<u32>().unwrap();
    let image_height: u32 = matches.value_of("HEIGHT")
        .unwrap().parse::<u32>().unwrap();
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
