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

fn image_correction(pixels: &[f32]) -> Vec<u8> {
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

fn supersampling(pixels: &[f32], w: usize, h: usize, ratio: usize) -> Vec<f32> {
    let mut result = vec![0.0; w * h * 4];
    for j in 0..h {
        for i in 0..w {
            let mut color = [0.0f32; 4];
            for k in 0..ratio*ratio {
                color[0] += pixels[4 * (w * ratio * (ratio * j + k / ratio) + ratio * i + k % ratio)];
                color[1] += pixels[4 * (w * ratio * (ratio * j + k / ratio) + ratio * i + k % ratio) + 1];
                color[2] += pixels[4 * (w * ratio * (ratio * j + k / ratio) + ratio * i + k % ratio) + 2];
                color[3] += pixels[4 * (w * ratio * (ratio * j + k / ratio) + ratio * i + k % ratio) + 3];
            }
            color[0] /= (ratio * ratio) as f32;
            color[1] /= (ratio * ratio) as f32;
            color[2] /= (ratio * ratio) as f32;
            color[3] /= (ratio * ratio) as f32;
            result[4 * (w*j + i)] = color[0];
            result[4 * (w*j + i) + 1] = color[1];
            result[4 * (w*j + i) + 2] = color[2];
            result[4 * (w*j + i) + 3] = color[3];
        }
    }
    result
}

use clap::{Arg, App};
use rand::prelude::*;

fn main() {
    let materials = vec![
        Material::complex(1.0, 0.5, 1.1,
                                Color::new(0.9, 0.32, 0.36),
                             Color::new(0.0, 0.0, 0.0),
                             0.0),
        Material::simple(Color::new(0.2, 0.76, 0.1),
                          Color::new(0.5, 0.5, 0.5),
                          32.0),
        Material::simple(Color::new(0.9, 0.9, 0.2),
                          Color::new(0.8, 0.8, 0.8),
                          32.0),
        Material::complex(1.0, 0.0, 1.1,
                          Color::new(1.0, 1.0, 1.0),
                         Color::new(0.0, 0.0, 0.0), 0.0),
        Material::with_texture(0, Color::zero(), 0.0),
        Material::complex(0.5, 0.0, 1.1,
                          Color::new(0.8, 0.8, 0.8),
                         Color::new(0.0, 0.0, 0.0),
                         0.0),
        Material::simple(Color::new(0.780392, 0.568627, 0.113725),
                          Color::new(0.992157, 0.941176, 0.807843),
                         27.8974),
        Material::with_texture_normal(1, 2, Color::zero(), 0.0),
        Material::with_texture_normal(3, 4, Color::zero(), 0.0),
        Material::with_emission(Color::new(1.0, 1.0, 1.0)),
        Material::with_texture(3, Color::zero(), 0.0),
    ];
    let textures = vec![
        "resources/checker.png".to_string(),
        "resources/brickwall.jpg".to_string(),
        "resources/brickwall_normal.jpg".to_string(),
        "resources/Brick_Wall_012_COLOR.jpg".to_string(),
        "resources/Brick_Wall_012_NORM.jpg".to_string()
    ];
    let mut scene = Scene::new(materials.clone(), textures.clone(),
        vec![
            Sphere {
                origin: Point3::new(-1.0, -2.0, -3.0),
                radius: 2.0,
                mat_id: 0
            },
            Sphere {
                origin: Point3::new(3.5, -2.0, -3.0),
                radius: 2.0,
                mat_id: 3
            },
            Sphere {
                origin: Point3::new(2.0, -2.0, -6.0),
                radius: 2.0,
                mat_id: 2
            },
            Sphere {
                origin: Point3::new(7.0, -2.0, -8.0),
                radius: 2.0,
                mat_id: 8
            },
        ],
        vec![],
        vec![
            PointLight {
                pos: Point3f::new(0.0, 5.0, 2.0),
                // pos: Point3f::new(3.0, -2.0, -1.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
        ],
        Point3::new(0.0, 2.0, 15.0)
    );
    let mut scene2 = Scene::new(materials.clone(), textures.clone(),
        vec![
            Sphere {
                origin: Point3::new(-1.0, -2.0, -3.0),
                radius: 2.0,
                mat_id: 0
            },
            Sphere {
                origin: Point3::new(3.5, -2.0, -3.0),
                radius: 2.0,
                mat_id: 3
            },
            Sphere {
                origin: Point3::new(2.0, -2.0, -6.0),
                radius: 2.0,
                mat_id: 2
            },
            Sphere {
                origin: Point3::new(7.0, -2.0, -8.0),
                radius: 2.0,
                mat_id: 8
            },
            Sphere {
                origin: Point3::new(0.0, 0.0, -11000.0),
                radius: 10000.0,
                mat_id: 9
            }
        ],
        vec![],
        vec![
            PointLight {
                pos: Point3f::new(0.0, 1000.0, 1000.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
            PointLight {
                pos: Point3f::new(0.0, -1000.0, 1000.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
        ],
        Point3::new(0.0, 2.0, 15.0)
    );
    let mut scene3 = Scene::new(materials.clone(), textures.clone(),
        vec![],
        vec![],
        vec![
            PointLight {
                pos: Point3f::new(0.0, 5.0, 2.0),
                // pos: Point3f::new(3.0, -2.0, -1.0),
                emission_color: Color::new(1.0, 1.0, 1.0)
            },
        ],
        Point3::new(0.0, 2.0, 15.0)
    );

    let mut scene5 = Scene::new(materials.clone(), textures.clone(),
       vec![
           Sphere {
               origin: Point3::new(-1.0, -2.0, -3.0),
               radius: 2.0,
               mat_id: 1
           },
           Sphere {
               origin: Point3::new(3.5, -2.0, -3.0),
               radius: 2.0,
               mat_id: 1
           },
           Sphere {
               origin: Point3::new(2.0, -2.0, -6.0),
               radius: 2.0,
               mat_id: 1
           },
           Sphere {
               origin: Point3::new(7.0, -2.0, -8.0),
               radius: 2.0,
               mat_id: 8
           },
       ],
       vec![],
       vec![
           PointLight {
               pos: Point3f::new(4.0, -1.5, -0.5),
               emission_color: Color::new(1.0, 1.0, 1.0)
           },
       ],
       Point3::new(0.0, 2.0, 15.0)
    );

    let matches = App::new("rusty-tracer")
        .version("1.0")
        .author("Philsik Chang <lasagnaphil@snu.ac.kr>")
        .about("A ray tracer written in rust")
        .arg(Arg::with_name("single-thread")
            .short("s").long("single-thread"))
        .arg(Arg::with_name("SCENE").index(1))
        .arg(Arg::with_name("WIDTH").index(2))
        .arg(Arg::with_name("HEIGHT").index(3))
        .arg(Arg::with_name("SAMPLE_RATE").index(4))
        .arg(Arg::with_name("OUTPUT").index(5))
        .get_matches();

    let scene_num: u32 = matches.value_of("SCENE")
        .unwrap_or("0").parse::<u32>().unwrap();
    let multithreading = !matches.is_present("single-thread");

    let scene = match scene_num {
        1 => &mut scene,
        2 => &mut scene2,
        3 | 4 => &mut scene3,
        5 => &mut scene5,
        _ => panic!("Invalid scene number!")
    };

    let filename = "resources/teapot.obj";
    let mut obj_meshes = obj_to_meshes(filename);
    println!("Model {} successfully loaded.", filename);
    for mut mesh in obj_meshes {
        let transform = Matrix4f::from_translation(
            Vector3f::new(-6.0, -3.0, -6.0));
        mesh.transform(transform);
        mesh.mat_id = 6;
        match scene_num {
            1 | 2 | 5 => { scene.add_mesh(mesh); },
            _ => {}
        }
    }

    let mut plane = Mesh::plane(4, 5.0);
    let transform = Matrix4f::from_translation(Vector3f::new(0.0, -4.0, 0.0));
    let transform = transform * Matrix4f::from_scale(100.0);
    plane.transform(transform);
    plane.create_tangents();
    match scene_num {
        2 => { plane.mat_id = 3; }
        _ => {}
    };
    scene.add_mesh(plane);

    let mut cube = Mesh::cube(4);
    let transform = Matrix4f::from_translation(Vector3f::new(0.0, 0.0, -15.0));
    let transform = transform * Matrix4f::from_scale(8.0);
    let transform = transform * Matrix4f::from_angle_y(cgmath::Deg(30.0));
    cube.transform(transform);
    cube.mat_id = 7;
    match scene_num {
        3 => { cube.mat_id = 8; },
        4 => { cube.mat_id = 10; }
        _ => {},
    };
    cube.create_tangents();
    scene.add_mesh(cube);

    let image_width: u32 = matches.value_of("WIDTH")
        .unwrap().parse::<u32>().unwrap();
    let image_height: u32 = matches.value_of("HEIGHT")
        .unwrap().parse::<u32>().unwrap();
    let sample_rate = matches.value_of("SAMPLE_RATE")
        .unwrap().parse::<u32>().unwrap();
    let output_filename = matches.value_of("OUTPUT")
        .unwrap();

    let trace_width = image_width * sample_rate;
    let trace_height = image_height * sample_rate;
    let mut pixels = vec![0.0; (trace_width * trace_height * 4) as usize];

    let time = Local::now();

    if multithreading {
        println!("Starting ray tracer with multiple threads.");
        let bands: Vec<(usize, &mut [f32])> = pixels
            .chunks_mut((trace_width * 4) as usize)
            .enumerate()
            .collect();

        let mut pb = Arc::new(Mutex::new(ProgressBar::new(bands.len() as u64)));
        pb.lock().unwrap().format("╢▌▌░╟");

        bands.into_par_iter().for_each(|(i, band)| {
            let pb = pb.clone();
            pb.lock().unwrap().inc();
            let band_upper_left = (0, i as u32);
            let band_lower_right = (trace_width, i as u32 + 1);
            scene.render(band, (trace_width, trace_height), band_upper_left, band_lower_right);
        });

        println!("Ray tracing complete!")
    }
    else {
        println!("Starting ray tracer with single thread.");

        let time = Local::now();
        scene.render(&mut pixels, (trace_width, trace_height), (0, 0), (trace_width, trace_height));

        println!("Ray tracing complete!")
    }

    println!("Elapsed time: {}ms", Local::now().signed_duration_since(time).num_milliseconds());

    if sample_rate > 1 {
        println!("Performing supersampling.");
        pixels = supersampling(&pixels, image_width as usize, image_height as usize, sample_rate as usize);
    }

    println!("Performing gamma correction.");
    let image_data = image_correction(&pixels);

    println!("Saving image.");
    image::save_buffer(&output_filename, &image_data, image_width, image_height, image::RGBA(8)).unwrap();

    println!("Rendered image saved at {}", &output_filename);
}
