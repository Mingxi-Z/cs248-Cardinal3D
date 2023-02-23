
#include "../lib/mathlib.h"
#include "debug.h"

bool BBox::hit(const Ray& ray, Vec2& times) const {

    // TODO (PathTracer):
    // Implement ray - bounding box intersection test
    // If the ray intersected the bounding box within the range given by
    // [times.x,times.y], update times with the new intersection times.
    Vec3 a = 1.0f / ray.dir;
    Vec3 b = -ray.point / ray.dir;

    Vec3 t_min_vec = a * min + b;
    Vec3 t_max_vec = a * max + b;

    for(int i = 0; i < 3; i++) {
        if(a[i] < 0) {
            std::swap(t_min_vec[i], t_max_vec[i]);
        }
    }

    float t_min = std::max({t_min_vec[0], t_min_vec[1], t_min_vec[2], ray.dist_bounds.x});
    float t_max = std::min({t_max_vec[0], t_max_vec[1], t_max_vec[2], ray.dist_bounds.y});
    
    if(t_min >= times.x && t_max <= times.y) {
        times.x = t_min;
        times.y = t_max;
    }

    if(t_max >= t_min) {
        return true;
    }

    return false;
}
