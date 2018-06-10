use super::shapes::*;
use super::scene::*;

use cgmath::prelude::*;

#[derive(Clone)]
pub struct BVH {
    pub bounding_volumes: Vec<BoundingVolume>
}

impl BVH {
    pub fn new(scene: &Scene) -> Self {
        let bounding_volumes: Vec<BoundingVolume> = scene.meshes.iter()
            .map(|m| BoundingVolume::from_mesh(m))
            .collect();
        BVH { bounding_volumes }
    }

    fn intersect(ray: &Ray, scene: &Scene) {
    }
}

pub struct Octree {
    children: Box<[Octree; 8]>
}

impl Octree {

}