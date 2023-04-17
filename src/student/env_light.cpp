
#include "../rays/env_light.h"
#include "debug.h"
#include "../lib/mathlib.h"

#include <limits>

namespace PT {

Light_Sample Env_Map::sample() const {

    Light_Sample ret;
    ret.distance = std::numeric_limits<float>::infinity();

    // TODO (PathTracer): Task 7
    // Uniformly sample the sphere. Tip: implement Samplers::Sphere::Uniform
    // Samplers::Sphere::Uniform uniform;
    // ret.direction = uniform.sample(ret.pdf);

    // Once you've implemented Samplers::Sphere::Image, remove the above and
    // uncomment this line to use importance sampling instead.
    ret.direction = sampler.sample(ret.pdf);
    ret.radiance = sample_direction(ret.direction);
    return ret;
}

Spectrum Env_Map::sample_direction(Vec3 dir) const {


    dir = -dir;
    dir.y = -dir.y;
    float theta = std::acos(dir.y);
    float phi = std::atan2(dir.z, dir.x);
    if (phi < 0.0f) {
        phi += 2.0f * PI_F;
    }

    float u = phi / (2.0f * PI_F);
    float v = 1.0f - theta / PI_F;
    auto [width, height] = image.dimension();

    int x = std::min(size_t(floor(u * width)), width - 1);
    int y = std::min(size_t(floor(v * height)), height - 1);

    //int x = static_cast<int>(u * width);
    //int y = static_cast<int>(v * height);
    int x1 = (x + 1) % width;
    int y1 = (y + 1) % height;

    float alpha = u * width - x;
    float beta = v * height - y;

    Spectrum s00 = image.at(x, y);
    Spectrum s01 = image.at(x, y1);
    Spectrum s10 = image.at(x1, y);
    Spectrum s11 = image.at(x1, y1);

    
    return s00 * (1 - alpha) * (1 - beta) + s10 * alpha * (1 - beta) + s01 * (1 - alpha) * beta + s11 * alpha * beta;
}

Light_Sample Env_Hemisphere::sample() const {
    Light_Sample ret;
    ret.direction = sampler.sample(ret.pdf);
    ret.radiance = radiance;
    ret.distance = std::numeric_limits<float>::infinity();
    return ret;
}

Spectrum Env_Hemisphere::sample_direction(Vec3 dir) const {
    if(dir.y > 0.0f) return radiance;
    return {};
}

Light_Sample Env_Sphere::sample() const {
    Light_Sample ret;
    ret.direction = sampler.sample(ret.pdf);
    ret.radiance = radiance;
    ret.distance = std::numeric_limits<float>::infinity();
    return ret;
}

Spectrum Env_Sphere::sample_direction(Vec3) const {
    return radiance;
}

} // namespace PT
