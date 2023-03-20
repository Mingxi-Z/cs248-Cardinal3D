
#include "../rays/bvh.h"
#include "debug.h"
#include <stack>

namespace PT {

// construct BVH hierarchy given a vector of prims
template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {

    // NOTE (PathTracer):
    // This BVH is parameterized on the type of the primitive it contains. This allows
    // us to build a BVH over any type that defines a certain interface. Specifically,
    // we use this to both build a BVH over triangles within each Tri_Mesh, and over
    // a variety of Objects (which might be Tri_Meshes, Spheres, etc.) in Pathtracer.
    //
    // The Primitive interface must implement these two functions:
    //      BBox bbox() const;
    //      Trace hit(const Ray& ray) const;
    // Hence, you may call bbox() and hit() on any value of type Primitive.

    // Keep these two lines of code in your solution. They clear the list of nodes and
    // initialize member variable 'primitives' as a vector of the scene prims
    nodes.clear();
    primitives = std::move(prims);

    // TODO (PathTracer): Task 3
    // Modify the code ahead to construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.
    //
    // Please use the SAH as described in class.  We recomment the binned build from lecture.
    // In general, here is a rough sketch:
    //
    //  For each axis X,Y,Z:
    //     Try possible splits along axis, evaluate SAH for each
    //  Take minimum cost across all axes.
    //  Partition primitives into a left and right child group
    //  Compute left and right child bboxes
    //  Make the left and right child nodes.
    //
    //
    // While a BVH is conceptually a tree structure, the BVH class uses a single vector (nodes)
    // to store all the nodes. Therefore, BVH nodes don't contain pointers to child nodes,
    // but rather the indices of the
    // child nodes in this array. Hence, to get the child of a node, you have to
    // look up the child index in this vector (e.g. nodes[node.l]). Similarly,
    // to create a new node, don't allocate one yourself - use BVH::new_node, which
    // returns the index of a newly added node.
    //
    // As an example of how to make nodes, the starter code below builds a BVH with a
    // root node that encloses all the primitives and its two descendants at Level 2.
    // For now, the split is hardcoded such that the first primitive is put in the left
    // child of the root, and all the other primitives are in the right child.
    // There are no further descendants.

    // edge case
    if(primitives.empty()) {
        return;
    }

    // compute bounding box for all primitives
    BBox bb;
    for(size_t i = 0; i < primitives.size(); ++i) {
        bb.enclose(primitives[i].bbox());
    }

    // set up root node (root BVH). Notice that it contains all primitives.
    size_t root_node_addr = new_node();
    Node& node = nodes[root_node_addr];
    node.bbox = bb;
    node.start = 0;
    node.size = primitives.size();

    // Creating Buckets in each direction
    const int n_buckets = 8;
    
    struct bucket {
        BBox bbox;
        size_t prim_count;
        bucket() : bbox(), prim_count(0) {}
    };

    // Create n_buckets buckets for each axis

    size_t iter = 0;

    do {
        if(nodes[iter].size <= max_leaf_size) {
            ++iter;
            continue;
        }
        
        // Get dimension of the box
        Vec3 min_pt = nodes[iter].bbox.min;
        Vec3 max_pt = nodes[iter].bbox.max;

        Vec3 dimension = max_pt - min_pt;

        int best_axis = 0;
        float best_split_coord = 0;

        float min_cost = std::numeric_limits<float>::max();

        // For each axis
        for(int axis = 0; axis < 3; ++axis) {
            bucket buckets[n_buckets];

            float bucket_width = dimension[axis] / n_buckets;

            // For each p in prims
            size_t end = nodes[iter].start + nodes[iter].size;

            for(size_t j = nodes[iter].start; j < end; ++j) {
                Primitive& p = primitives[j];

                // Compute Bucket Index for p
                BBox pbox = p.bbox();
                Vec3 center = pbox.center();

                int bucket_idx = (int)std::floorf((center[axis] - min_pt[axis]) / bucket_width);

                bucket_idx = std::clamp(bucket_idx, 0, n_buckets - 1);

                // Add to bucket
                buckets[bucket_idx].bbox.enclose(pbox);
                buckets[bucket_idx].prim_count++;
            }

            // Evaluate SAH for each plane
            bucket left{}, right{};
            float parent_surface = nodes[iter].bbox.surface_area();

            for(int l = 0; l < n_buckets - 1; ++l) {
                bucket& cur = buckets[l];
                // Skip empty bucket
                if(cur.prim_count == 0) continue;

                left.bbox.enclose(cur.bbox);
                left.prim_count += cur.prim_count;

                right.bbox.reset();
                right.prim_count = 0;

                for(int r = l + 1; r < n_buckets; ++r) {
                    cur = buckets[r];
                    right.bbox.enclose(cur.bbox);
                    right.prim_count += cur.prim_count;
                }

                float cost = left.bbox.surface_area() / parent_surface * left.prim_count;
                cost += right.bbox.surface_area() / parent_surface * right.prim_count;

                if(cost < min_cost) {
                    best_axis = axis;
                    min_cost = cost;
                    best_split_coord = bucket_width * (l + 1);
                }
            }
        }

        best_split_coord += min_pt[best_axis];

        // Reorder nodes based on split_coordinate and axis
        size_t right_start = std::partition(primitives.begin() + nodes[iter].start,
                                            primitives.begin() + nodes[iter].start + nodes[iter].size,
                                            [best_axis, best_split_coord](const Primitive& x) {
                                                BBox cur = x.bbox();
                                                return cur.center()[best_axis] <= best_split_coord;
                                            }) -
                             primitives.begin();

        if(right_start == nodes[iter].start || right_start == nodes[iter].start + nodes[iter].size) {
            ++iter;
            continue;
        }

        // Construct Left and Right node
        BBox left{}, right{};
        int left_cnt = 0, right_cnt = 0;
        for(size_t i = nodes[iter].start; i < right_start; ++i) {
            left.enclose(primitives[i].bbox());
            left_cnt++;
        }

        for(size_t i = right_start; i < nodes[iter].start + nodes[iter].size; ++i) {
            right.enclose(primitives[i].bbox());
            right_cnt++;
        }

        assert(left.min.x >= nodes[iter].bbox.min.x && left.max.y <= nodes[iter].bbox.max.y &&
               right.min.x >= nodes[iter].bbox.min.x && right.max.y <= nodes[iter].bbox.max.y);
        nodes[iter].l = new_node(left, nodes[iter].start, left_cnt, 0, 0);
        nodes[iter].r = new_node(right, right_start, right_cnt, 0, 0);
        ++iter;
    } while(iter < nodes.size());
}

template<typename Primitive>
Trace BVH<Primitive>::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

    Trace ret;
    
    auto find_hit = [&, this, ray](auto& self, size_t idx, auto& nodes, Trace& ret) -> void { 
        Node cur = nodes[idx];

        if(nodes[idx].is_leaf()) {
            size_t start = cur.start;
            size_t end = start + cur.size;
            for(size_t i = cur.start; i < end; ++i) {
                Trace hit = primitives[i].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            Vec2 time1, time2;
            bool hit1 = nodes[cur.l].bbox.hit(ray, time1);
            bool hit2 = nodes[cur.r].bbox.hit(ray, time2);

            size_t first = cur.l,
                second = cur.r;

            if(hit1 && hit2) {
                if(time1.x >= time2.x) {
                    std::swap(first, second);
                    std::swap(time1, time2);
                }
                self(self, first, nodes, ret);

                if(time2.x < ray.dist_bounds.y) {
                    self(self, second, nodes, ret);
                }
            } else if(hit1){
                self(self, first, nodes, ret);
            } else if(hit2) {
                self(self, second, nodes, ret);
            }
        }
    };

    Vec2 time;
    if (nodes.empty())
        return ret;
    bool hit = nodes[0].bbox.hit(ray, time);
    if(!hit) return ret;

    find_hit(find_hit, 0, nodes, ret);
    
    return ret;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    build(std::move(prims), max_leaf_size);
}

template<typename Primitive>
BVH<Primitive> BVH<Primitive>::copy() const {
    BVH<Primitive> ret;
    ret.nodes = nodes;
    ret.primitives = primitives;
    ret.root_idx = root_idx;
    return ret;
}

template<typename Primitive>
bool BVH<Primitive>::Node::is_leaf() const {
    return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
    Node n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.l = l;
    n.r = r;
    nodes.push_back(n);
    return nodes.size() - 1;
}

template<typename Primitive>
BBox BVH<Primitive>::bbox() const {
    return nodes[root_idx].bbox;
}

template<typename Primitive>
std::vector<Primitive> BVH<Primitive>::destructure() {
    nodes.clear();
    return std::move(primitives);
}

template<typename Primitive>
void BVH<Primitive>::clear() {
    nodes.clear();
    primitives.clear();
}

template<typename Primitive>
size_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                                 const Mat4& trans) const {

    std::stack<std::pair<size_t, size_t>> tstack;
    tstack.push({root_idx, 0});
    size_t max_level = 0;

    if(nodes.empty()) return max_level;

    while(!tstack.empty()) {

        auto [idx, lvl] = tstack.top();
        max_level = std::max(max_level, lvl);
        const Node& node = nodes[idx];
        tstack.pop();

        Vec3 color = lvl == level ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(1.0f);
        GL::Lines& add = lvl == level ? active : lines;

        BBox box = node.bbox;
        box.transform(trans);
        Vec3 min = box.min, max = box.max;

        auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

        edge(min, Vec3{max.x, min.y, min.z});
        edge(min, Vec3{min.x, max.y, min.z});
        edge(min, Vec3{min.x, min.y, max.z});
        edge(max, Vec3{min.x, max.y, max.z});
        edge(max, Vec3{max.x, min.y, max.z});
        edge(max, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

        if(node.l && node.r) {
            tstack.push({node.l, lvl + 1});
            tstack.push({node.r, lvl + 1});
        } else {
            for(size_t i = node.start; i < node.start + node.size; i++) {
                size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
                max_level = std::max(c, max_level);
            }
        }
    }
    return max_level;
}

} // namespace PT
