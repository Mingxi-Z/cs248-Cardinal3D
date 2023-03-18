
#include <assimp/Importer.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../util/camera.h"
#include "../util/rand.h"
#include "../rays/samplers.h"
#include "debug.h"

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

    float sensor_x = (-half_width) * (1 - screen_coord.x) + half_width * screen_coord.x;
    float sensor_y = (-half_height) * (1 - screen_coord.y) + half_height * screen_coord.y;


    Vec3 sensor_pos(sensor_x, sensor_y, -1.0f);
    // Tip: compute the ray direction in view space and use
    // the camera space to world space transform (iview) to transform the ray back into world space.
    Ray r(Vec3(0, 0, 0), sensor_pos);

    if (aperture > 0) {
        // Sample point on lens
        auto sampleDisk = [](float radius) -> Vec2 {
            float angle = RNG::unit() * 2 * PI_F;
            float r = radius * RNG::unit();
            return Vec2(r * cosf(angle), r * sinf(angle));
        };
        Vec2 pLens = sampleDisk(aperture);

        // Compute point on plane of focus
        float ft = focal_dist / r.dir.z;
        Vec3 pFocus = r.at(ft);

        // Update ray for effect of lens
        r.point = Vec3(pLens.x, pLens.y, 0.0f);
        r.dir = (r.point - pFocus);
    }

    r.transform(iview);
    return r;
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
        element.curvature *= 0.001;
        element.thickness *= 0.001;
        element.eta *= 0.001;
        element.aperture = element.aperture * 0.001 / 2.0;
        lens_elements.push_back(element);
    }

    input_file.close();
}
