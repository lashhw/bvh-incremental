#include <iostream>
#include <bvh/triangle.hpp>
#include <bvh/bvh.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/primitive_intersectors.hpp>
#include "happly.h"

typedef bvh::Bvh<float> bvh_t;
typedef bvh::Triangle<float> triangle_t;
typedef bvh::Vector3<float> vector_t;
typedef bvh::Ray<float> ray_t;
typedef bvh::Bvh<float>::Node node_t;
typedef bvh::SweepSahBuilder<bvh_t> builder_t;
typedef bvh::SingleRayTraverser<bvh_t> traverser_t;
typedef bvh::ClosestPrimitiveIntersector<bvh_t, triangle_t> intersector_t;
typedef traverser_t::Statistics statistics_t;

constexpr int num_bits = 6;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "usage: ./bvh_incremental MODEL RAY" << std::endl;
        exit(EXIT_FAILURE);
    }

    happly::PLYData ply_data(argv[1]);
    std::vector<std::array<double, 3>> v_pos = ply_data.getVertexPositions();
    std::vector<std::vector<size_t>> f_idx = ply_data.getFaceIndices<size_t>();

    std::vector<triangle_t> triangles;
    for (auto &face : f_idx) {
        triangles.emplace_back(vector_t(v_pos[face[0]][0], v_pos[face[0]][1], v_pos[face[0]][2]),
                               vector_t(v_pos[face[1]][0], v_pos[face[1]][1], v_pos[face[1]][2]),
                               vector_t(v_pos[face[2]][0], v_pos[face[2]][1], v_pos[face[2]][2]));
    }

    auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(triangles.data(), triangles.size());
    auto global_bbox = bvh::compute_bounding_boxes_union(bboxes.get(), triangles.size());
    std::cout << "global_bbox = ("
              << global_bbox.min[0] << ", " << global_bbox.min[1] << ", " << global_bbox.min[2] << "), ("
              << global_bbox.max[0] << ", " << global_bbox.max[1] << ", " << global_bbox.max[2] << ")" << std::endl;

    std::cout << "building..." << std::endl;
    bvh_t full_bvh;
    builder_t builder(full_bvh);
    builder.build(global_bbox, bboxes.get(), centers.get(), triangles.size());
    std::cout << "node_count = " << full_bvh.node_count << std::endl;

    std::cout << "copying..." << std::endl;
    bvh_t reduced_bvh;
    reduced_bvh.node_count = full_bvh.node_count;
    reduced_bvh.nodes = std::make_unique<node_t[]>(reduced_bvh.node_count);
    std::copy(full_bvh.nodes.get(),
              full_bvh.nodes.get() + full_bvh.node_count,
              reduced_bvh.nodes.get());
    reduced_bvh.primitive_indices = std::make_unique<size_t[]>(triangles.size());
    std::copy(full_bvh.primitive_indices.get(),
              full_bvh.primitive_indices.get() + triangles.size(),
              reduced_bvh.primitive_indices.get());

    std::cout << "updating..." << std::endl;
    std::stack<size_t> stack;
    stack.push(0);
    while (!stack.empty()) {
        const auto &node = reduced_bvh.nodes[stack.top()];
        stack.pop();
        if (node.is_leaf())
            continue;
        for (int i = 0; i < 2; i++) {
            size_t child_idx = node.first_child_or_primitive + i;
            auto &child_node = reduced_bvh.nodes[child_idx];
            stack.push(child_idx);
            for (int j = 0; j < 3; j++) {
                float min = node.bounds[2 * j];
                float max = node.bounds[2 * j + 1];
                float diff = max - min;
                int exp;
                frexpf(diff, &exp);
                int k = exp + 1;
                float &child_min = child_node.bounds[2 * j];
                float &child_max = child_node.bounds[2 * j + 1];
                assert(child_min >= min);
                assert(child_max <= max);
                int r = (int)((child_min - min) / ldexpf(1.0f, k) * ldexpf(1.0f, num_bits));
                int s = (int)((child_max - min) / ldexpf(1.0f, k) * ldexpf(1.0f, num_bits));
                float new_min = min + ldexpf((float)r, k - num_bits);
                float new_max = min + ldexpf((float)(s + 1), k - num_bits);
                assert(new_min <= child_min);
                assert(new_max >= child_max);
                child_min = new_min;
                child_max = new_max;
            }
        }
    }

    std::cout << "traversing..." << std::endl;
    traverser_t full_traverser(full_bvh);
    traverser_t reduced_traverser(reduced_bvh);
    intersector_t full_intersector(full_bvh, triangles.data());
    intersector_t reduced_intersector(reduced_bvh, triangles.data());
    std::ifstream ray_file(argv[2], std::ios::in | std::ios::binary);
    float r[7];
    int count = 0;
    int same_count = 0;
    std::ofstream full_trv_file("full_trv.bin", std::ios::out | std::ios::binary);
    std::ofstream full_ist_file("full_ist.bin", std::ios::out | std::ios::binary);
    std::ofstream reduced_trv_file("reduced_trv.bin", std::ios::out | std::ios::binary);
    std::ofstream reduced_ist_file("reduced_ist.bin", std::ios::out | std::ios::binary);
    for (; ray_file.read((char*)(&r), 7 * sizeof(float)); count++) {
        ray_t ray(
            vector_t(r[0], r[1], r[2]),
            vector_t(r[3], r[4], r[5]),
            0.f,
            r[6]
        );
        statistics_t full_stat;
        statistics_t reduced_stat;
        auto full_result = full_traverser.traverse(ray, full_intersector, full_stat);
        auto reduced_result = reduced_traverser.traverse(ray, reduced_intersector, reduced_stat);
        bool same = true;
        if (full_result) {
            if (reduced_result) {
                same &= (full_result->primitive_index == reduced_result->primitive_index);
                same &= (full_result->intersection.t == reduced_result->intersection.t);
                same &= (full_result->intersection.u == reduced_result->intersection.u);
                same &= (full_result->intersection.v == reduced_result->intersection.v);
            } else {
                same = false;
            }
        } else {
            if (reduced_result) {
                same = false;
            }
        }
        if (same)
            same_count++;
        full_trv_file.write((char*)(&full_stat.traversal_steps), sizeof(full_stat.traversal_steps));
        full_ist_file.write((char*)(&full_stat.intersections), sizeof(full_stat.intersections));
        reduced_trv_file.write((char*)(&reduced_stat.traversal_steps), sizeof(reduced_stat.traversal_steps));
        reduced_ist_file.write((char*)(&reduced_stat.intersections), sizeof(reduced_stat.intersections));
    }

    std::cout << same_count << "/" << count << std::endl;
}
