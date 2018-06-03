use super::shapes::*;
use super::scene::*;

struct BVH {
    bounding_volumes: Vec<BoundingVolume>
}

impl BVH {
    pub fn new(meshes: &[Mesh]) -> Self {
        let bounding_volumes: Vec<BoundingVolume> = meshes.iter()
            .map(|m| m.get_bounding_volume())
            .collect();
        BVH { bounding_volumes }
    }

    fn intersect(ray: &Ray) -> Color {

    }
}

pub struct Octree {
    children: Box<[Octree; 8]>
}

impl Octree {

}