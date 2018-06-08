extern crate cgmath;

use super::bvh::*;

use std::mem;
use std::f32;

use cgmath::prelude::*;
use cgmath::Point2;
use cgmath::Point3;
use cgmath::Vector2;
use cgmath::Vector3;
use cgmath::Matrix4;

pub type Point2f = Point2<f32>;
pub type Point3f = Point3<f32>;
pub type Vector2f = Vector2<f32>;
pub type Vector3f = Vector3<f32>;
pub type Matrix4f = Matrix4<f32>;
pub type Color = Vector3<f32>;

pub struct Ray {
    pub origin: Point3f,
    pub dir: Vector3f
}

impl Ray {
    pub fn new(origin: Point3f, dir: Vector3f) -> Self {
        Ray { origin: origin, dir: dir }
    }
}

pub trait Intersectable {
    type IntersectResult;
    fn intersect(self: &Self, ray: &Ray) -> Option<Self::IntersectResult>;
}

#[derive(Clone)]
pub struct Sphere {
    pub origin: Point3f,
    pub radius: f32,
    pub mat_id: usize
}

impl Intersectable for Sphere {
    type IntersectResult = f32;
    fn intersect(self: &Self, ray: &Ray) -> Option<f32> {
        let disp = self.origin - ray.origin;
        let ip = ray.dir.dot(disp);
        if ip < 0.0 { return None; }
        let discriminant = ip * ip - disp.magnitude2() + self.radius * self.radius;
        if discriminant >= 0.0 {
            let tnear = ip - discriminant.sqrt();
            let tfar = ip + discriminant.sqrt();
            if tnear > 1e-4 { Some(tnear) } else { Some(tfar) }
        }
        else {
            None
        }
    }
}

#[derive(Clone)]
pub struct Vertex {
    pub pos: Point3f,
    pub normal: Vector3f,
    pub tex: Point2f
}

impl Vertex {
    pub fn empty() -> Vertex {
        Vertex::from_floats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }
    pub fn new(pos: Point3f, normal: Vector3f, tex: Point2f) -> Vertex {
        Vertex { pos, normal, tex }
    }
    pub fn from_floats(x: f32, y: f32, z: f32, nx: f32, ny: f32, nz: f32, u: f32, v: f32) -> Vertex {
        Vertex { pos: Point3f::new(x, y, z), normal: Vector3f::new(nx, ny, nz), tex: Point2f::new(u, v) }
    }
}

#[derive(Clone)]
pub struct Triangle {
    pub vertices: [Vertex; 3],
}

impl Triangle {
    pub fn from_vertices(vertices: &[Vertex]) -> Self {
        Triangle { vertices: [vertices[0].clone(), vertices[1].clone(), vertices[2].clone()] }
    }

    pub fn normal_of_pos(self: &Self, pos: Point3f) -> Vector3f { let e0 = self.vertices[1].pos - self.vertices[0].pos;
        let e1 = self.vertices[2].pos - self.vertices[1].pos;
        let e2 = self.vertices[0].pos - self.vertices[2].pos;
        let a = e0.cross(pos - self.vertices[1].pos).magnitude();
        let b = e1.cross(pos - self.vertices[2].pos).magnitude();
        let c = e2.cross(pos - self.vertices[0].pos).magnitude();
        (a * self.vertices[0].normal + b * self.vertices[0].normal + c * self.vertices[0].normal) / (a + b + c)
    }
}

impl Intersectable for Triangle {
    type IntersectResult = (f32, f32, f32);
    // Using Moller-Trumbore intersection algorithm
    // Note: Quite dirty!
    fn intersect(self: &Self, ray: &Ray) -> Option<(f32, f32, f32)> {
        const EPSILON: f32 = 0.0000001;
        let edge1 = self.vertices[1].pos - self.vertices[0].pos;
        let edge2 = self.vertices[2].pos - self.vertices[0].pos;
        let h = ray.dir.cross(edge2);
        let a = edge1.dot(h);
        if a > -EPSILON && a < EPSILON { return None; }
        let f = 1.0/a;
        let s = ray.origin - self.vertices[0].pos;
        let u = f * s.dot(h);
        if u < 0.0 || u > 1.0 { return None; }
        let q = s.cross(edge1);
        let v = f * ray.dir.dot(q);
        if v < 0.0 || u + v > 1.0 { return None; }
        let t = f * edge2.dot(q);
        if t > EPSILON { Some((t, u, v)) }
            else { None }
    }
}

#[derive(Clone)]
pub struct Mesh {
    pub vertices: Vec<Vertex>,
    pub mat_id: usize
}

impl Mesh {
    pub fn from_vertices(vertices: Vec<Vertex>, mat_id: usize) -> Self {
        Mesh { vertices, mat_id }
    }

    pub fn cube(mat_id: usize) -> Mesh {
        Mesh::from_vertices(
            vec![
                Vertex::from_floats(-0.5, -0.5, -0.5,  0.0,  0.0, -1.0,  0.0, 0.0),
                Vertex::from_floats(0.5, -0.5, -0.5,  0.0,  0.0, -1.0,  1.0, 0.0),
                Vertex::from_floats(0.5,  0.5, -0.5,  0.0,  0.0, -1.0,  1.0, 1.0),
                Vertex::from_floats(0.5,  0.5, -0.5,  0.0,  0.0, -1.0,  1.0, 1.0),
                Vertex::from_floats(-0.5,  0.5, -0.5,  0.0,  0.0, -1.0,  0.0, 1.0),
                Vertex::from_floats(-0.5, -0.5, -0.5,  0.0,  0.0, -1.0,  0.0, 0.0),
                Vertex::from_floats(-0.5, -0.5,  0.5,  0.0,  0.0, 1.0,   0.0, 0.0),
                Vertex::from_floats(0.5, -0.5,  0.5,  0.0,  0.0, 1.0,   1.0, 0.0),
                Vertex::from_floats(0.5,  0.5,  0.5,  0.0,  0.0, 1.0,   1.0, 1.0),
                Vertex::from_floats(0.5,  0.5,  0.5,  0.0,  0.0, 1.0,   1.0, 1.0),
                Vertex::from_floats(-0.5,  0.5,  0.5,  0.0,  0.0, 1.0,   0.0, 1.0),
                Vertex::from_floats(-0.5, -0.5,  0.5,  0.0,  0.0, 1.0,   0.0, 0.0),
                Vertex::from_floats(-0.5,  0.5,  0.5, -1.0,  0.0,  0.0,  1.0, 0.0),
                Vertex::from_floats(-0.5,  0.5, -0.5, -1.0,  0.0,  0.0,  1.0, 1.0),
                Vertex::from_floats(-0.5, -0.5, -0.5, -1.0,  0.0,  0.0,  0.0, 1.0),
                Vertex::from_floats(-0.5, -0.5, -0.5, -1.0,  0.0,  0.0,  0.0, 1.0),
                Vertex::from_floats(-0.5, -0.5,  0.5, -1.0,  0.0,  0.0,  0.0, 0.0),
                Vertex::from_floats(-0.5,  0.5,  0.5, -1.0,  0.0,  0.0,  1.0, 0.0),
                Vertex::from_floats(0.5,  0.5,  0.5,  1.0,  0.0,  0.0,  1.0, 0.0),
                Vertex::from_floats(0.5,  0.5, -0.5,  1.0,  0.0,  0.0,  1.0, 1.0),
                Vertex::from_floats(0.5, -0.5, -0.5,  1.0,  0.0,  0.0,  0.0, 1.0),
                Vertex::from_floats(0.5, -0.5, -0.5,  1.0,  0.0,  0.0,  0.0, 1.0),
                Vertex::from_floats(0.5, -0.5,  0.5,  1.0,  0.0,  0.0,  0.0, 0.0),
                Vertex::from_floats(0.5,  0.5,  0.5,  1.0,  0.0,  0.0,  1.0, 0.0),
                Vertex::from_floats(-0.5, -0.5, -0.5,  0.0, -1.0,  0.0,  0.0, 1.0),
                Vertex::from_floats(0.5, -0.5, -0.5,  0.0, -1.0,  0.0,  1.0, 1.0),
                Vertex::from_floats(0.5, -0.5,  0.5,  0.0, -1.0,  0.0,  1.0, 0.0),
                Vertex::from_floats(0.5, -0.5,  0.5,  0.0, -1.0,  0.0,  1.0, 0.0),
                Vertex::from_floats(-0.5, -0.5,  0.5,  0.0, -1.0,  0.0,  0.0, 0.0),
                Vertex::from_floats(-0.5, -0.5, -0.5,  0.0, -1.0,  0.0,  0.0, 1.0),
                Vertex::from_floats(-0.5,  0.5, -0.5,  0.0,  1.0,  0.0,  0.0, 1.0),
                Vertex::from_floats(0.5,  0.5, -0.5,  0.0,  1.0,  0.0,  1.0, 1.0),
                Vertex::from_floats(0.5,  0.5,  0.5,  0.0,  1.0,  0.0,  1.0, 0.0),
                Vertex::from_floats(0.5,  0.5,  0.5,  0.0,  1.0,  0.0,  1.0, 0.0),
                Vertex::from_floats(-0.5,  0.5,  0.5,  0.0,  1.0,  0.0,  0.0, 0.0),
                Vertex::from_floats(-0.5, 0.5, -0.5, 0.0, 1.0, 0.0, 0.0, 1.0),
            ],
            mat_id
        )
    }
    pub fn plane(mat_id: usize) -> Mesh {
        Mesh::from_vertices(
            vec![
                Vertex::from_floats(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
                Vertex::from_floats(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0),
                Vertex::from_floats(1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0),
                Vertex::from_floats(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0),
            ],
            mat_id
        )
    }

    pub fn transform(self: &mut Self, transform: Matrix4f) {
        for v in &mut self.vertices {
            v.pos = transform.transform_point(v.pos);
            v.normal = transform.transform_vector(v.normal);
        }
    }

    pub fn transform_with(self: &mut Self, transform_fn: fn(&mut Vertex) -> ()) {
        for vertex in &mut self.vertices {
            transform_fn(vertex);
        }
    }
}


pub const NUM_PLANE_NORMALS: usize = 9;

lazy_static! {
    pub static ref PLANE_NORMALS: [Vector3f; NUM_PLANE_NORMALS] = [
        Vector3f::new(1.0, 0.0, 0.0),
        Vector3f::new(0.0, 1.0, 0.0),
        Vector3f::new(0.0, 0.0, 1.0),
        Vector3f::new(1.0, 1.0, 0.0).normalize(),
        Vector3f::new(1.0, 0.0, 1.0).normalize(),
        Vector3f::new(0.0, 1.0, 1.0).normalize(),
        Vector3f::new(1.0, -1.0, 0.0).normalize(),
        Vector3f::new(1.0, 0.0, -1.0).normalize(),
        Vector3f::new(0.0, 1.0, -1.0).normalize(),
    ];
}

pub struct BoundingVolume {
    pub ranges: [(f32, f32); NUM_PLANE_NORMALS],
    pub mesh: Mesh
}

impl BoundingVolume {
    pub fn from_mesh(mesh: &Mesh) -> Self {
        let mut ranges = [(f32::INFINITY, -f32::INFINITY); NUM_PLANE_NORMALS];
        for i in 0..NUM_PLANE_NORMALS {
            let plane_normal = PLANE_NORMALS[i];
            for vertex in &mesh.vertices {
                let dist = plane_normal.dot(vertex.pos - Point3f::origin());
                if dist < ranges[i].0 { ranges[i].0 = dist; }
                if dist > ranges[i].1 { ranges[i].1 = dist; }
            }
        }
        BoundingVolume { ranges, mesh: mesh.clone() }
    }

    pub fn intersect(self: &Self, ray: &Ray, normal_orig_angles: &[f32], normal_dir_angles: &[f32])
        -> Option<(f32, f32, usize)> {
        let mut tnear_final = -f32::INFINITY;
        let mut tfar_final = f32::INFINITY;
        let mut imin = NUM_PLANE_NORMALS;
        for (i, range) in self.ranges.iter().enumerate() {
            let tnear = (range.1 - normal_orig_angles[i]) / normal_dir_angles[i];
            let tfar = (range.0 - normal_orig_angles[i]) / normal_dir_angles[i];
            let (tnear, tfar) = if normal_dir_angles[i] < 0.0 { (tnear, tfar) } else { (tfar, tnear) };
            if tnear > tnear_final {
                tnear_final = tnear;
                imin = i;
            }
            if tfar < tfar_final {
                tfar_final = tfar;
            }
            if tnear_final > tfar_final { return None; }
        }
        if imin < NUM_PLANE_NORMALS {
            Some((tnear_final, tfar_final, imin))
                /*
            if tnear_final > 0.0 {
                Some((tnear_final, tfar_final, imin))
            }
            else if tfar_final > 0.0 {
                Some((0.0, tfar_final, imin))
            }
            else {
                None
            }
            */
        } else { None }
    }
}

pub enum Shape {
    Sphere(Sphere),
    Triangle(Triangle),
    BoundingVolume(usize)
}

