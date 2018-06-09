
use std::path::Path;
use tobj;
use super::shapes::*;
use cgmath::prelude::*;

pub fn obj_to_meshes(filename: &str) -> Vec<Mesh> {
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
        let has_indices = mesh.indices.len() > 0;

        println!("has_normals: {}", has_normals);
        println!("has_texcoords: {}", has_texcoords);
        println!("has_indices: {}", has_indices);

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

        if has_indices {
            if !has_normals {
                for i in 0..(num_indices / 3) {
                    let e1 = vertices[mesh.indices[3 * i + 1] as usize].pos - vertices[mesh.indices[3 * i] as usize].pos;
                    let e2 = vertices[mesh.indices[3 * i + 2] as usize].pos - vertices[mesh.indices[3 * i] as usize].pos;
                    let normal = e1.cross(e2).normalize();
                    vertices[mesh.indices[3 * i] as usize].normal += normal;
                    vertices[mesh.indices[3 * i + 1] as usize].normal += normal;
                    vertices[mesh.indices[3 * i + 2] as usize].normal += normal;
                }
                for i in 0..num_vertices {
                    vertices[i].normal = vertices[i].normal.normalize();
                }
            }

            vertices = mesh.indices.iter()
                .map(|ind| vertices[*ind as usize].clone())
                .collect();
        }

        Mesh { vertices, tangents: None, mat_id }
    }).collect();

    meshes
}