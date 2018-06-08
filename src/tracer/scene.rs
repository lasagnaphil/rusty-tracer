use tobj;
use super::shapes::*;
use super::bvh::*;

use std;
use std::f32;
use std::collections::HashMap;

use cgmath::prelude::*;
use cgmath::Point3;
use cgmath::Vector2;
use cgmath::Vector3;

#[derive(Clone)]
pub struct Material {
    pub surface_color: Color,
    pub reflectivity: f32,
    pub transparency: f32,
    pub refractive_index: f32,
    pub specular_color: Color,
    pub shininess: f32
}

pub struct PointLight {
    pub pos: Point3f,
    pub emission_color: Color
}

pub struct Scene {
    pub camera_pos: Point3f,
    pub materials: Vec<Material>,
    pub spheres: Vec<Sphere>,
    pub meshes: Vec<Mesh>,
    pub point_lights: Vec<PointLight>,
    pub bvh: Option<BVH>
}

const BIAS: f32 = 1e-4;

impl Scene {
    const MAX_RAY_DEPTH: u32 = 5;

    pub fn new(materials: Vec<Material>, spheres: Vec<Sphere>, meshes: Vec<Mesh>,
               point_lights: Vec<PointLight>, camera_pos: Point3f) -> Self {
        Scene { materials, spheres, meshes, point_lights, camera_pos, bvh: None }
    }


    pub fn add_mesh(self: &mut Self, mesh: Mesh) {
        self.meshes.push(mesh);
    }

    fn get_material(self: &Self, id: usize) -> &Material {
        &self.materials[id]
    }

    fn get_material_mut(self: &mut Self, id: usize) -> &mut Material {
        &mut self.materials[id]
    }

    fn get_triangles<'a>(self: &'a Self) -> Vec<Triangle> {
        self.meshes.iter()
            .flat_map(|m| {
                m.vertices
                    .chunks(3)
                    .map(|c| Triangle::from_vertices(c))
            }).collect()
    }

    pub fn build_bvh(self: &mut Self) {
        let bvh = BVH::new(self);
        self.bvh = Some(bvh);
    }

    fn trace(self: &Self, ray: &Ray, depth: u32) -> Color {
        let mut closest_shape: Option<Shape> = None;
        let mut tnear = f32::INFINITY;
        let mut u = 0.0;
        let mut v = 0.0;
        let mut mat_id = 0;

        // Find the nearest collision of ray with scene
        for sphere in &self.spheres {
            if let Some(t) = sphere.intersect(&ray) {
                if t < (tnear - BIAS) {
                    tnear = t;
                    closest_shape = Some(Shape::Sphere(sphere.clone()));
                    mat_id = sphere.mat_id;
                }
            }
        }

        if let Some(bvh) = self.bvh.as_ref() {
            let mut normal_orig_angles = [0.0; NUM_PLANE_NORMALS];
            let mut normal_dir_angles = [0.0; NUM_PLANE_NORMALS];
            for i in 0..NUM_PLANE_NORMALS {
                normal_orig_angles[i] = PLANE_NORMALS[i].dot(ray.origin - Point3f::origin());
                normal_dir_angles[i] = PLANE_NORMALS[i].dot(ray.dir);
                // scene.meshes
            }

            for (i, bounding_volume) in bvh.bounding_volumes.iter().enumerate() {
                if let Some((tn, tf, plane_i)) = bounding_volume.intersect(
                    ray, &normal_orig_angles, &normal_dir_angles) {
                    if tn < (tnear - BIAS) {
                        tnear = tn;
                        closest_shape = Some(Shape::BoundingVolume(i));
                        mat_id = bounding_volume.mesh.mat_id;
                    }
                }
            }

            if let Some(Shape::BoundingVolume(bv_index)) = closest_shape {
                tnear = f32::INFINITY;
                let mesh = &bvh.bounding_volumes[bv_index].mesh;
                for i in 0..(mesh.vertices.len() / 3) {
                    let triangle = Triangle::from_vertices(&mesh.vertices[3 * i..3 * i + 3]);
                    if let Some((t, tu, tv)) = triangle.intersect(&ray) {
                        if t < (tnear - BIAS) {
                            tnear = t;
                            u = tu;
                            v = tv;
                            closest_shape = Some(Shape::Triangle(triangle));
                            mat_id = mesh.mat_id;
                        }
                    }
                }
                // TODO:
                // Right now this is an ugly hack
                // It doesn't catch some edge cases, so there will be more shadow rays
                // This will be solved when using BV tree instead of BV list
                if let Some(Shape::Triangle(_)) = closest_shape {} else { return Color::zero(); }
            }
        }
        else {
            for mesh in &self.meshes {
                for i in 0..(mesh.vertices.len() / 3) {
                    let triangle = Triangle::from_vertices(&mesh.vertices[3 * i..3 * i + 3]);
                    if let Some((t, tu, tv)) = triangle.intersect(&ray) {
                        if t < (tnear - BIAS) {
                            tnear = t;
                            u = tu;
                            v = tv;
                            closest_shape = Some(Shape::Triangle(triangle));
                            mat_id = mesh.mat_id;
                        }
                    }
                }
            }
        }

        if let Some(shape) = closest_shape {
            self.trace_shape(shape, mat_id, tnear, u, v, ray, depth)
        } else {
            Color::zero()
        }
    }

    fn trace_shape(self: &Self, closest_shape: Shape,
                   mat_id: usize, tnear: f32, u: f32, v: f32, ray: &Ray, depth: u32) -> Color {
        let mat = self.get_material(mat_id);
        let hit_pos = ray.origin + tnear * ray.dir;
        let hit_normal = match closest_shape {
            Shape::Sphere(sphere) => {
                (hit_pos - sphere.origin).normalize()
            },
            Shape::Triangle(triangle) => {
                let normal = (1.0 - u - v) * triangle.vertices[0].normal
                    + u * triangle.vertices[1].normal
                    + v * triangle.vertices[2].normal;
                reflect(ray.dir, normal)
            },
            Shape::BoundingVolume(bounding_volume) => {
                eprintln!("Unreachable code pos!");
                std::process::exit(2);
            }
        };
        let incident_angle = -ray.dir.dot(hit_normal);
        let hit_normal = if incident_angle > 0.0 { hit_normal } else { -hit_normal };
        if (mat.transparency > 0.0 || mat.reflectivity > 0.0) && depth < Scene::MAX_RAY_DEPTH {
            let n = if incident_angle > 0.0 {
                1.0 / mat.refractive_index
            } else {
                mat.refractive_index
            };
            let r0 = ((n - 1.0) / (n + 1.0)).powi(2);
            let fresnel = r0 + (1.0 - r0) * (1.0 - incident_angle.abs()).powi(5);

            let reflect_dir = reflect(ray.dir, hit_normal);
            let reflection_ray = Ray::new(hit_pos + hit_normal * BIAS, reflect_dir);
            let reflection_color = self.trace(&reflection_ray, depth + 1);

            let refraction_color = if mat.transparency > 0.0 {
                let refract_dir = refract(ray.dir, hit_normal, n);
                let refraction_ray = Ray::new(hit_pos - hit_normal * BIAS, refract_dir);
                self.trace(&refraction_ray, depth + 1)
            } else { Color::zero() };

            mat.surface_color.mul_element_wise(
                reflection_color * fresnel + refraction_color * (1.0 - fresnel) * mat.transparency
            )
        } else {
            let mut surface_color = Color::zero();
            for point_light in &self.point_lights {
                if point_light.emission_color != Color::zero() {
                    let shadow_ray = Ray::new(hit_pos + hit_normal * BIAS, (point_light.pos - hit_pos).normalize());
                    let mut is_shadow = false;
                    for other_sphere in &self.spheres {
                        if other_sphere.intersect(&shadow_ray).is_some() {
                            is_shadow = true; break;
                        }
                    }
                    for mesh in &self.meshes {
                        for i in 0..(mesh.vertices.len()/3) {
                            let triangle = Triangle::from_vertices(&mesh.vertices[3*i..3*i+3]);
                            if triangle.intersect(&shadow_ray).is_some() {
                                is_shadow = true; break;
                            }
                        }
                    }
                    if !is_shadow {
                        let shadow_angle = hit_normal.dot(shadow_ray.dir);
                        surface_color += shadow_angle.max(0.0) *
                            mat.surface_color.mul_element_wise(point_light.emission_color);
                        let specular_angle = -ray.dir.dot(shadow_ray.dir);
                        surface_color += specular_angle.max(0.0).powf(mat.shininess) *
                            mat.specular_color.mul_element_wise(point_light.emission_color);
                    };
                }
            }
            surface_color
        }
    }

    pub fn render(self: &Self, pixels: &mut [f32], bounds: (u32, u32), upper_left: (u32, u32), lower_right: (u32, u32)) {
        let fov = 60.0f32;
        let tangent = (fov / 2.0f32).to_radians().tan();
        let aspect_ratio = (bounds.1-1) as f32 / (bounds.0-1) as f32;

        let iter_width = lower_right.0 - upper_left.0;
        let iter_height = lower_right.1 - upper_left.1;
        for j in 0..iter_height {
            for i in 0..iter_width {
                let prim_ray_dir = Vector3::new(
                    ((2.0 * (upper_left.0 + i) as f32 / (bounds.0 - 1) as f32) - 1.0) * tangent,
                    -((2.0 * (upper_left.1 + j) as f32 / (bounds.0 - 1) as f32) - aspect_ratio) * tangent,
                    -1.0).normalize();

                let prim_ray = Ray::new(self.camera_pos, prim_ray_dir);

                let color = self.trace(&prim_ray, 0);
                pixels[(4*(j * iter_width + i)) as usize] = color.x;
                pixels[(4*(j * iter_width + i) + 1) as usize] = color.y;
                pixels[(4*(j * iter_width + i) + 2) as usize] = color.z;
                pixels[(4*(j * iter_width + i) + 3) as usize] = 1.0;
            }
        }
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
