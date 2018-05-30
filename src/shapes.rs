extern crate cgmath;

use cgmath::prelude::*;
use cgmath::Point3;
use cgmath::Vector2;
use cgmath::Vector3;

pub type Point3f = Point3<f32>;
pub type Vector2f = Vector2<f32>;
pub type Vector3f = Vector3<f32>;

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
    type IntersectResult = (f32, f32);
    fn intersect(self: &Self, ray: &Ray) -> Option<(f32, f32)> {
        let disp = self.origin - ray.origin;
        let ip = ray.dir.dot(disp);
        if ip < 0.0 { return None; }
        let discriminant = ip * ip - disp.magnitude2() + self.radius * self.radius;
        if discriminant >= 0.0 {
            Some((ip - discriminant.sqrt(), ip + discriminant.sqrt()))
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
    pub tex: Vector2f
}

pub struct Triangle {
    pub vertices: [Vertex; 3],
    pub mat_id: usize
}

impl Traceable for Triangle {
    type IntersectResult = f32;
    // Using Moller-Trumbore intersection algorithm
    // Note: Quite dirty!
    fn intersect(self: &Self, ray: &Ray) -> Option<f32> {
        const EPSILON: f32 = 0.0000001;
        let edge1 = self.vertices[1].pos - self.vertices[0].pos;
        let edge2 = self.vertices[2].pos - self.vertices[0].pos;
        let h = ray.dir.cross(edge2);
        let a = edge1.dot(h);
        if a > -EPSILON && a < EPSILON { return None; }
        let f = 1.0/a;
        let s = ray.origin - self.vertices[1].pos;
        let u = f * s.dot(h);
        if u < 0.0 || u > 1.0 { return None; }
        let q = s.cross(edge1);
        let v = f * ray.dir.dot(q);
        if v < 0.0 || u + v > 1.0 { return None; }
        let t = f * edge2.dot(q);
        if t > EPSILON { Some(t) }
            else { None }
    }
}