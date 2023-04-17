
#include "../rays/samplers.h"
#include "../util/rand.h"
#include "debug.h"

namespace Samplers {

Vec2 Rect::Uniform::sample(float& pdf) const {

    // TODO (PathTracer): Task 1
    // Generate a uniformly random point on a rectangle of size (size.x by size.y)

    // Tip: consider using RNG::unit()  (see util/rand.h)
    Vec2 xy(size.x * RNG::unit(), size.y * RNG::unit());

    // PDF is the probability density of the chosen sample
    // the PDF should integrate to 1 over the whole rectangle
    pdf = 1.0f / (size.x * size.y); 

    // Return the randomly generated point
    return xy;
}

Vec3 Hemisphere::Cosine::sample(float& pdf) const {

    // TODO (PathTracer): Task 6
    // You may implement this, but don't have to.
    return Vec3();
}

Vec3 Sphere::Uniform::sample(float& pdf) const {

    // TODO (PathTracer): Task 7
    // Generate a uniformly random point on the unit sphere (or equivalently, direction)
    // Tip: start with Hemisphere::Uniform
    Hemisphere::Uniform sampler;

    Vec3 dir = sampler.sample(pdf);

    // Randomly flip the direction
    if (RNG::coin_flip()) {
        dir = -dir;
    }

    // Set the PDF to be equal for all directions
    pdf *= 0.5f / PI_F;

    return dir;
}

Sphere::Image::Image(const HDR_Image& image) {

    // TODO (PathTracer): Task 7
    // Set up importance sampling for a spherical environment map image.

    // You may make use of the pdf, cdf, and total members, or create your own
    // representation.

    const auto [_w, _h] = image.dimension();
    w = _w;
    h = _h;
    float total = 0.0f;

    // Compute the pdf and cdf arrays for the image
    pdf.reserve(w * h);
    cdf.reserve(w * h);
    for (int j = 0; j < (int) h; ++j) {
        for (int i = 0; i < (int) w; ++i) {
            Spectrum color = image.at(i, j);
            float luminance = color.luma();
            pdf.push_back(luminance);
            total += luminance;
            if (i == 0 && j == 0) {
                cdf.push_back(luminance);
            } else {
                cdf.push_back(cdf.back() + luminance);
            }
        }
    }

    // Normalize pdf and cdf arrays
    const float inv_total = 1.0f / total;
    for (auto& p : pdf) {
        p *= inv_total;
    }
    for (auto& p : cdf) {
        p *= inv_total;
    }
}

Vec3 Sphere::Image::sample(float& out_pdf) const {

    // TODO (PathTracer): Task 7
    // Use your importance sampling data structure to generate a sample direction.
    // Tip: std::upper_bound can easily binary search your CDF

    out_pdf = 1.0f; // what was the PDF (again, PMF here) of your chosen sample?

    float r = RNG::unit();
    auto it = std::upper_bound(cdf.begin(), cdf.end(), r);
    int index = it - cdf.begin();

    int i = index % w;
    int j = index / w;
    float u = (i + 0.5f) / w;
    float v = (j + 0.5f) / h;
    float phi = 2 * PI_F * u;
    float cos_theta = 2 * v - 1;
    float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
    Vec3 dir(sin_theta * cos(phi), 2 * v - 1, sin_theta * sin(phi));

    out_pdf = pdf[index] * w * h / (2 * PI_F * PI_F * sin_theta);

    return dir;
}

Vec3 Point::sample(float& pmf) const {

    pmf = 1.0f;
    return point;
}

Vec3 Two_Points::sample(float& pmf) const {
    if(RNG::unit() < prob) {
        pmf = prob;
        return p1;
    }
    pmf = 1.0f - prob;
    return p2;
}

Vec3 Hemisphere::Uniform::sample(float& pdf) const {

    float Xi1 = RNG::unit();
    float Xi2 = RNG::unit();

    float theta = std::acos(Xi1);
    float phi = 2.0f * PI_F * Xi2;

    float xs = std::sin(theta) * std::cos(phi);
    float ys = std::cos(theta);
    float zs = std::sin(theta) * std::sin(phi);

    pdf = 1.0f / (2.0f * PI_F);
    return Vec3(xs, ys, zs);
}

} // namespace Samplers
