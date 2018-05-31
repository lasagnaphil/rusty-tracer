extern crate cgmath;

use std::mem;

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

pub struct Ray {
    pub origin: Point3f,
    pub dir: Vector3f
}

impl Ray {
    pub fn new(origin: Point3f, dir: Vector3f) -> Self {
        Ray { origin: origin, dir: dir }
    }
}

pub trait Traceable {
    type IntersectResult;
    fn intersect(self: &Self, ray: &Ray) -> Option<Self::IntersectResult>;
}

#[derive(Clone)]
pub struct Sphere {
    pub origin: Point3f,
    pub radius: f32,
    pub mat_id: usize
}

impl Traceable for Sphere {
    type IntersectResult = f32;
    fn intersect(self: &Self, ray: &Ray) -> Option<f32> {
        let disp = self.origin - ray.origin;
        let ip = ray.dir.dot(disp);
        if ip < 0.0 { return None; }
        let discriminant = ip * ip - disp.magnitude2() + self.radius * self.radius;
        if discriminant >= 0.0 {
            let tnear = ip - discriminant.sqrt();
            let tfar = ip + discriminant.sqrt();
            if tnear < 0.0 { Some(tfar) } else { Some(tnear) }
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
    fn empty() -> Vertex {
        Vertex::from_floats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }
    fn new(pos: Point3f, normal: Vector3f, tex: Point2f) -> Vertex {
        Vertex { pos, normal, tex }
    }
    fn from_floats(x: f32, y: f32, z: f32, nx: f32, ny: f32, nz: f32, u: f32, v: f32) -> Vertex {
        Vertex { pos: Point3f::new(x, y, z), normal: Vector3f::new(nx, ny, nz), tex: Point2f::new(u, v) }
    }
}

#[derive(Clone)]
pub struct Triangle {
    pub vertices: [Vertex; 3],
    pub mat_id: usize
}

impl Triangle {
    pub fn from_vertices(vertices: &[Vertex], mat_id: usize) -> Self {
        unsafe {
            let mut triangle_vertices: [Vertex; 3] = mem::uninitialized();
            triangle_vertices.clone_from_slice(vertices);
            Triangle { vertices: triangle_vertices, mat_id }
        }
    }

    pub fn normal_of_pos(self: &Self, pos: Point3f) -> Vector3f {
        let e0 = self.vertices[1].pos - self.vertices[0].pos;
        let e1 = self.vertices[2].pos - self.vertices[1].pos;
        let e2 = self.vertices[0].pos - self.vertices[2].pos;
        let a = e0.cross(pos - self.vertices[1].pos).magnitude();
        let b = e1.cross(pos - self.vertices[2].pos).magnitude();
        let c = e2.cross(pos - self.vertices[0].pos).magnitude();
        (a * self.vertices[0].normal + b * self.vertices[0].normal + c * self.vertices[0].normal) / (a + b + c)
    }

    pub fn transform(self: Self, transform: Matrix4f) -> Self {
        let mut vertices = self.vertices;
        for v in &mut vertices {
            v.pos = transform.transform_point(v.pos);
            v.normal = transform.transform_vector(v.normal);
        }
        Triangle {
            vertices, mat_id: self.mat_id
        }
    }
}

impl Traceable for Triangle {
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

pub struct Model {
    pub triangles: Vec<Triangle>,
    pub mat_id: usize
}

impl Model {
    pub fn from_vertices(vertices: Vec<Vertex>, mat_id: usize) -> Self {
        let triangles: Vec<Triangle> = vertices
            .chunks(3)
            .map(|c| Triangle::from_vertices(c, mat_id))
            .collect();
        Model { triangles, mat_id }
    }

    pub fn cube(mat_id: usize) -> Model {
        Model::from_vertices(
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
}

pub enum Shape {
    Sphere(Sphere),
    Triangle(Triangle),
}

impl Shape {
    pub fn mat_id(self: &Self) -> usize {
        match self {
            Shape::Sphere(sphere) => sphere.mat_id,
            Shape::Triangle(triangle) => triangle.mat_id
        }
    }
}