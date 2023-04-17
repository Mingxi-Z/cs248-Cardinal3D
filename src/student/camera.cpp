
#include <assimp/Importer.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../util/camera.h"
#include "../util/rand.h"
#include "../rays/samplers.h"
#include "../rays/shapes.h"
#include "debug.h"
#include "../rays/bsdf.h"
#include "../lib/bbox.h"
#include "../util/thread_pool.h"

#define LENS_DIST 0.018f
#define SENSOR_SIZE 0.135f / 3
#define EXIT_PUPIL_SAMPLE 100

PT::Trace spherical_element_hit (float radius, const Ray &ray);
Vec3 refract (Vec3 &out_dir, Vec3 &normal, float eta, bool &was_internal);
void get_cardinal_point (const Ray& in, const Ray& out, float &pz, float&fz);

Ray Camera::generate_ray(Vec2 screen_coord) const {

    // TODO (PathTracer): Task 1
    //
    // The input screen_coord is a normalized screen coordinate [0,1]^2
    //
    // You need to transform this 2D point into a 3D position on the sensor plane, which is
    // located one unit away from the pinhole in camera space (aka view space).
    //
    // You'll need to compute this position based on the vertial field of view
    // (vert_fov) of the camera, and the aspect ratio of the output image (aspect_ratio).
    Vec2 top_left, top_right, bottom_left, bottom_right;

    float half_fov = Radians(get_fov() / 2);

    float half_height = tanf(half_fov);
    float half_width = half_height * get_ar();

    // Calculate the film plane for 135mm camera
    float aspect = half_height / half_width;
    float diagonal = SENSOR_SIZE;
    float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
    float y = aspect * x;
    
    half_height = y / 2;
    half_width = x / 2;

    float sensor_x = (-half_width) * (1 - screen_coord.x) + half_width * screen_coord.x;
    float sensor_y = (-half_height) * (1 - screen_coord.y) + half_height * screen_coord.y;


    Vec3 sensor_pos(-sensor_x, -sensor_y, 0.0f);
    // Tip: compute the ray direction in view space and use
    // the camera space to world space transform (iview) to transform the ray back into world space.
    Ray r(Vec3(0, 0, 0), sensor_pos.unit());

    if (lens_elements.empty()){
        r.transform(iview);
        return r;
    }
    
    Lens_Element last = lens_elements.at(lens_elements.size() - 1);

    auto sampleDisk = [](float radius) -> Vec2 {
        float angle = RNG::unit() * 2 * PI_F;
        float r = radius * sqrt(RNG::unit());
        return Vec2(r * cosf(angle), r * sinf(angle));
    };
    Vec2 pLens2 = sampleDisk(last.aperture);

    Vec3 dest(pLens2.x, pLens2.y, -last.thickness);

    Vec3 orig(-sensor_x, -sensor_y, 0.0f);
    
    r = Ray(orig, (dest - orig).unit()); 
    
    bool success = false;
    Ray camera = trace_lens_ray(r, &success);
    //Ray camera = trace_lens_ray(r, &success);
    //camera.dir = -camera.dir;
    camera.transform(iview);
    // camera.dir = -camera.dir;
    if (!success)
        camera.dist_bounds.y = EPS_F;

    return camera;
}

// Real Camera
void Camera::load_lens(std::string lens_path) {
    std::cout << lens_path << std::endl;
    std::ifstream input_file(lens_path);
    std::string line;

    if (!input_file.is_open()) {
        // Error opening file
        return;
    }

    if (!lens_elements.empty())
        lens_elements.clear();

    while (std::getline(input_file, line)) {
        if (line[0] == '#') {
            continue; // Skip lines that start with #
        }
        std::istringstream iss(line);
        Lens_Element element;
        if (!(iss >> element.curvature >> element.thickness >> element.eta >> element.aperture)) {
            std::cout << "Error: could not parse line \"" << line << "\"" << std::endl;
            break;
        }
        element.curvature *= 0.001f;
        element.thickness *= 0.001f;
        element.aperture = element.aperture * 0.001f / 2.0f;
        lens_elements.emplace_back(element);
    }

    lens_elements.back().thickness = thick_lens_focus();
    input_file.close();
}

float Camera::lens_rear_z() const {
    assert (!lens_elements.empty());
    return lens_elements.back().thickness;
}

float Camera::lens_front_z() const {
    assert (!lens_elements.empty());
    float front_z = 0.0f;
    for (const auto& element : lens_elements) {
        front_z += element.thickness;
    }
    return front_z;
}

float Camera::rear_radius() const {
    return lens_elements.back().aperture;
}

Ray Camera::trace_lens_ray (const Ray &camera_ray, bool *success) const {
    float element_z = 0.0f;
    Ray lens_ray = camera_ray;
    PT::Trace result;
    Lens_Element elem;

    for (int i = lens_elements.size() - 1; i >= 0; i--) {         
        elem = lens_elements.at(i);
        element_z -= elem.thickness;
        
        float t = 0.0f;
        // Compute intersection of ray with lens element
        if(elem.curvature == 0.0f) {
            // Case when lens element is flat (aperture stop)
            float a = 1.0f / lens_ray.dir.z;
            float b = -lens_ray.point.z / lens_ray.dir.z;
            t = a * element_z + b;
        } else {
            // Case when lens element is spherical
            float radius = elem.curvature;
            float z_center = element_z + elem.curvature;
            
            // Shift ray origin to center of sphere for intersection test
            lens_ray.point.z -= z_center;
            result = spherical_element_hit (radius, lens_ray);
            if (!result.hit) {
                *success = false;
                return lens_ray;
            }
            lens_ray.point.z += z_center;
            t = result.distance;
        }
        // Test intersection point against element aperture
        Vec3 hit_location = lens_ray.at(t);
        float r2 = hit_location.x * hit_location.x + hit_location.y * hit_location.y;
        float ap = (elem.curvature == 0.0f) ? aperture : elem.aperture;
        if (r2 > ap * ap) {
            *success = false;
            return lens_ray;
        }
        lens_ray.point = hit_location;
        // Update ray path for element interface interaction
        if (elem.curvature != 0.0f) {
            // Case when lens element is spherical
            Vec3 w;
            float etaI = elem.eta;
            float etaT = (i > 0 && lens_elements[i - 1].eta != 0.0f) ?
                        lens_elements[i - 1].eta : 1.0f;
            bool was_internal = false;

            float eta = etaI / etaT;
            Vec3 temp = -lens_ray.dir;

            // Flip surface normal if ray is coming from inside lens element   
            if (dot(temp, result.normal) < 0.0f)
                result.normal = -result.normal;
            // Compute refracted ray direction
            Vec3 ref_dir = refract(temp, result.normal, eta, was_internal);
            if (was_internal){
                *success = false;
                return lens_ray;
            }
            lens_ray.dir = ref_dir;
        }
    }
    *success = true;
    return lens_ray;
}

Ray Camera::trace_lens_ray_reverse (const Ray &camera_ray, bool *success) const {
    float element_z = -lens_front_z ();
    Ray lens_ray = camera_ray;
    PT::Trace result;
    Lens_Element elem;

    for (int i = 0; i < (int) lens_elements.size(); i++) {         
        elem = lens_elements.at(i);
        float t = 0.0f;
        // Compute intersection of ray with lens element
        if(elem.curvature == 0.0f) {
            float a = 1.0f / lens_ray.dir.z;
            float b = -lens_ray.point.z / lens_ray.dir.z;
            t = a * element_z + b;
        } else {
            float radius = elem.curvature;
            float z_center = element_z + elem.curvature;
            
            lens_ray.point.z -= z_center;
            result = spherical_element_hit (radius, lens_ray);
            if (!result.hit) {
                *success = false;
                return lens_ray;
            }
            
            lens_ray.point.z += z_center;
            t = result.distance;
        }
        // Test intersection point against element aperture
        Vec3 hit_location = lens_ray.at(t);
        float r2 = hit_location.x * hit_location.x + hit_location.y * hit_location.y;
        if (r2 > elem.aperture * elem.aperture) {
            *success = false;
            return lens_ray;
        }
        lens_ray.point = hit_location;
        // Update ray path for element interface interaction
        if (elem.curvature == 0.0f) {
            element_z += elem.thickness;
            continue;
        }
            
        Vec3 w;
        float etaI = (i == 0 || lens_elements[i - 1].eta == 0)  ?
            1.0f : lens_elements[i - 1].eta;
        float etaT = (elem.eta != 0.0f) ? elem.eta : 1.0f;
        bool was_internal = false;

        float eta = etaI / etaT;
        Vec3 temp = -lens_ray.dir;
        if (dot(temp, result.normal) < 0.0f)
            result.normal = -result.normal;
        Vec3 ref_dir = refract(temp, result.normal, eta, was_internal);
        if (was_internal){
            *success = false;
            return lens_ray;
        }
        lens_ray.dir = ref_dir;
        element_z += elem.thickness;
    }
    *success = true;
    return lens_ray;
}

void Camera::thick_lens_approx (Vec2 &pz, Vec2 &fz) const {
    float x = .001f * SENSOR_SIZE;
    Ray in(Vec3(x, 0.0f, -lens_front_z() - 1.0f), Vec3(0.0f, 0.0f, 1.0f));
    bool success = false;
    Ray out = trace_lens_ray_reverse(in, &success);

    if (!success)
        return;

    get_cardinal_point(in, out, pz[0], fz[0]);

    out = Ray(Vec3(x, 0.0f, lens_rear_z() + 1.0f), Vec3(0, 0, -1));
    in = trace_lens_ray(out, &success);
    
    if (!success)
        return;

    get_cardinal_point(out, in, pz[1], fz[1]);
}

float Camera::thick_lens_focus (void) {
    Vec2 pz, fz;
    thick_lens_approx(pz, fz);

    float f = fz[0] - pz[0];
    float z = -focal_dist;
    float delta = 0.5f * (pz[1] - z + pz[0] -
        std::sqrt((pz[1] - z - pz[0]) * (pz[1] - z - 4 * f - pz[0])));
    return lens_rear_z() + delta;
}

PT::Trace spherical_element_hit (float radius, const Ray &ray) {
    PT::Trace ret;
    ret.origin = ray.point;
    ret.hit = false;       // was there an intersection?
    ret.distance = 0.0f;   // at what distance did the intersection occur?
    ret.position = Vec3{}; // where was the intersection?
    ret.normal = Vec3{};   // what was the surface normal at the intersection?

    float a = dot(-ray.point, ray.dir);
    float b = a * a - ray.point.norm_squared() + radius * radius;
    
    if(b < 0.0f) {
        return ret;
    }
    float t1 = a - std::sqrt(b);
    float t2 = a + std::sqrt(b);

    if(t1 > ray.dist_bounds.y || t2 < ray.dist_bounds.x) {
        return ret;
    }

    bool use_closer_t = (ray.dir.z > 0) ^ (radius < 0);
    float t = use_closer_t ? t1 : t2;
    ret.hit = true;
    ret.position = ray.at(t);
    ret.normal = (ret.position).unit();
    ret.distance = t;
    return ret;
}

Vec3 refract (Vec3 &out_dir, Vec3 &normal, float eta, bool &was_internal) {

    Vec3 in_dir;
    float cos_theta_out = dot(normal, out_dir);
    float sin2_theta_out = std::clamp(1.f - cos_theta_out * cos_theta_out, 0.0f, 1.0f);
    float sin2_theta_in = eta * eta * sin2_theta_out;

    if (sin2_theta_in >= 1) {
        was_internal = true;
        return in_dir;
    }

    float cos_theta_in = std::sqrt(1 - sin2_theta_in);

    in_dir = eta * -out_dir + (eta * cos_theta_out - cos_theta_in) * normal;
    was_internal = false;

    return in_dir;
}

void get_cardinal_point (const Ray& in, const Ray& out, float &pz, float&fz) {
    float a = in.point.x / out.dir.x;
    float b = out.point.x / out.dir.x;

    float tf = -b;
    fz = out.at(tf).z;
    float tp = a - b;
    pz = out.at(tp).z;
}