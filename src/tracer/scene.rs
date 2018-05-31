use tobj;
use super::shapes::*;

use std::f32;
use std::path::Path;

use cgmath::prelude::*;
use cgmath::Point3;
use cgmath::Vector2;
use cgmath::Vector3;

pub type Color = Vector3<f32>;

#[derive(Clone)]
pub struct Material {
    pub surface_color: Color,
    pub reflectivity: f32,
    pub transparency: f32,
    pub refractive_index: f32,
    pub emission_color: Color
}

pub struct Scene {
    pub camera_pos: Point3f,
    pub materials: Vec<Material>,
    pub spheres: Vec<Sphere>,
    pub meshes: Vec<Mesh>,
}

impl Scene {
    const MAX_RAY_DEPTH: u32 = 5;

    pub fn from_obj(filename: &str) -> Self {
        let obj = tobj::load_obj(&Path::new(filename));
        let (models, materials) = obj.unwrap();
        let meshes: Vec<Mesh> = models.iter().map(|model| {
            let mesh = &model.mesh;
            let mat_id = mesh.material_id.unwrap_or(0);
            let num_vertices = mesh.positions.len() / 3;
            let num_indices = mesh.indices.len();
            let mut vertices: Vec<Vertex> = vec![];

            let has_normals = mesh.normals.len() > 0;
            let has_texcoords = mesh.texcoords.len() > 0;

            for i in 0..num_vertices {
                let position = Point3f::new(mesh.positions[3*i], mesh.positions[3*i+1], mesh.positions[3*i+2]);
                let normal = if has_normals {
                    Vector3f::new(mesh.normals[3*i], mesh.normals[3*i+1], mesh.normals[3*i+2])
                } else {
                    Vector3f::zero()
                };
                let texcoords = if has_texcoords {
                    Point2f::new(mesh.texcoords[2*i], mesh.texcoords[2*i+1])
                } else {
                    Point2f::new(0.0, 0.0)
                };
                vertices.push(Vertex::new(position, normal, texcoords));
            }

            if mesh.indices.len() > 0 {
                vertices = mesh.indices.iter()
                    .map(|ind| vertices[*ind as usize].clone())
                    .collect();
            }

            let mut triangles: Vec<Triangle> = vertices
                .chunks(3)
                .map(|t| {
                    let mut vertices = [Vertex::empty(), Vertex::empty(), Vertex::empty()];
                    vertices.clone_from_slice(t);
                    Triangle { vertices, mat_id }
                })
                .collect();

            if !has_normals {
                for i in 0..(num_indices / 3) {
                    let e1 = triangles[i].vertices[1].pos - triangles[i].vertices[0].pos;
                    let e2 = triangles[i].vertices[2].pos - triangles[i].vertices[0].pos;
                    triangles[i].vertices[0].normal = e2.cross(e1).normalize();
                }
            }

            Mesh { triangles, mat_id }
        }).collect();

        let materials = vec![
            Material {
                surface_color: Color::new(0.8, 0.8, 0.8),
                reflectivity: 0.7,
                transparency: 0.0,
                refractive_index: 1.1,
                emission_color: Color::new(0.0, 0.0, 0.0)
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
                emission_color: Color::new(0.5, 0.5, 0.5)
            }
        ];

        Scene {
            meshes,
            materials,
            spheres: vec![
                Sphere {
                    origin: Point3::new(0.0, 0.0, -10020.0),
                    radius: 10000.0,
                    mat_id: 1
                },
                Sphere {
                    origin: Point3::new(0.0, 0.0, 10020.0),
                    radius: 10000.0,
                    mat_id: 1
                },
                Sphere {
                    origin: Point3::new(0.0, -10004.0, 0.0),
                    radius: 10000.0,
                    mat_id: 2
                }
            ],
            camera_pos: Point3f::new(0.0, 1.0, 10.0)
        }
    }

    fn get_material(self: &Self, id: usize) -> &Material {
        &self.materials[id]
    }

    fn get_material_mut(self: &mut Self, id: usize) -> &mut Material {
        &mut self.materials[id]
    }

    fn get_triangles(self: &Self) -> impl Iterator<Item = &Triangle> {
        self.meshes.iter().flat_map(|m| &m.triangles)
    }

    fn trace(self: &Self, ray: &Ray, depth: u32) -> Color {
        let mut closest_shape: Option<Shape> = None;
        let mut tnear = f32::INFINITY;
        let mut u = 0.0;
        let mut v = 0.0;

        // Find the nearest collision of ray with scene
        for sphere in &self.spheres {
            if let Some(t) = sphere.intersect(&ray) {
                if t < tnear {
                    tnear = t;
                    closest_shape = Some(Shape::Sphere(sphere.clone()));
                }
            }
        }
        for mesh in &self.meshes {
            for triangle in &mesh.triangles {
                if let Some((t, tu, tv)) = triangle.intersect(&ray) {
                    if t < tnear {
                        tnear = t; u = tu; v = tv;
                        closest_shape = Some(Shape::Triangle(triangle.clone()));
                    }
                }
            }
        }

        let shape = if let Some(shape) = closest_shape {
            shape
        } else {
            return Color::zero();
        };

        let mat = self.get_material(shape.mat_id());
        let hit_pos = ray.origin + tnear * ray.dir;
        let hit_normal = match shape {
            Shape::Sphere(sphere) => {
                (hit_pos - sphere.origin).normalize()
            },
            Shape::Triangle(triangle) => {
                let normal = u * triangle.vertices[0].normal
                    + v * triangle.vertices[1].normal
                    + (1.0 - u - v) * triangle.vertices[2].normal;
                reflect(ray.dir, normal)
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
            let reflection_ray = Ray::new(hit_pos + hit_normal * 1e-4, reflect_dir);
            let reflection_color = self.trace(&reflection_ray, depth + 1);

            let refraction_color = if mat.transparency > 0.0 {
                let refract_dir = refract(ray.dir, hit_normal, n);
                let refraction_ray = Ray::new(hit_pos - hit_normal * 1e-4, refract_dir);
                self.trace(&refraction_ray, depth + 1)
            } else { Color::zero() };

            mat.emission_color + mat.surface_color.mul_element_wise(
                reflection_color * fresnel + refraction_color * (1.0 - fresnel) * mat.transparency
            )
        } else {
            let mut surface_color = Color::zero();
            for light_sphere in &self.spheres {
                let light_mat = self.get_material(light_sphere.mat_id);
                if light_mat.emission_color != Color::zero() {
                    let shadow_ray = Ray::new(hit_pos, (light_sphere.origin - hit_pos).normalize());
                    let is_shadow = self.spheres.iter().any(|other_sphere| {
                        other_sphere.intersect(&shadow_ray).is_some()
                    });
                    let is_shadow = is_shadow || self.get_triangles().any(|triangle| {
                        triangle.intersect(&shadow_ray).is_some()
                    });
                    if !is_shadow {
                        let shadow_angle = hit_normal.dot(shadow_ray.dir);
                        if shadow_angle > 0.0 {
                            surface_color += shadow_angle *
                                mat.surface_color.mul_element_wise(light_mat.emission_color);
                        }
                    };
                }
            }

            mat.emission_color + surface_color
        }
    }

    pub fn render(self: &Self, pixels: &mut [f32], bounds: (u32, u32), upper_left: (u32, u32), lower_right: (u32, u32)) {
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
