
#include "../rays/shapes.h"
#include "debug.h"

namespace PT {

const char* Shape_Type_Names[(int)Shape_Type::count] = {"None", "Sphere"};

BBox Sphere::bbox() const {

    BBox box;
    box.enclose(Vec3(-radius));
    box.enclose(Vec3(radius));
    return box;
}

Trace Sphere::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!

    Trace ret;
    ret.origin = ray.point;
    ret.hit = false;       // was there an intersection?
    ret.distance = 0.0f;   // at what distance did the intersection occur?
    ret.position = Vec3{}; // where was the intersection?
    ret.normal = Vec3{};   // what was the surface normal at the intersection?

    float a = dot(-ray.point, ray.dir);
    float b = a * a - ray.point.norm_squared() + radius * radius;
    float t1 = a - std::sqrt(b);
    float t2 = a + std::sqrt(b);
    
    if(b < 0.0f) {
        return ret;
    }

    if(t1 > ray.dist_bounds.y || t2 < ray.dist_bounds.x) {
        return ret;
    }

    float t = t1;

    if(t1 < ray.dist_bounds.x) {
        if(t2 > ray.dist_bounds.x && t2 < ray.dist_bounds.y) {
            t = t2;
        } else {
            return ret;
        }
    }
    //ray.dist_bounds.y = t;

    ret.hit = true;
    ret.position = ray.at(t);
    ret.distance = t;
    ret.normal = (ret.position - bbox().center()).unit();
    return ret;
}

} // namespace PT
