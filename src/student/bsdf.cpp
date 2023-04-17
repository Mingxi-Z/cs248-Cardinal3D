
#include "../rays/bsdf.h"
#include "../util/rand.h"
#include "debug.h"

namespace PT {

Vec3 reflect(Vec3 dir) {

    // TODO (PathTracer): Task 6
    // Return reflection of dir about the surface normal (0,1,0).
    Vec3 normal(0.0f, 1.0f, 0.0f);
    Vec3 out = - dir + 2 * (dir * normal) * normal;
    return out;
}

Vec3 refract(Vec3 out_dir, float index_of_refraction, bool& was_internal) {

    // TODO (PathTracer): Task 6
    // Use Snell's Law to refract out_dir through the surface
    // Return the refracted direction. Set was_internal to false if
    // refraction does not occur due to total internal reflection,
    // and true otherwise.

    // When dot(out_dir,normal=(0,1,0)) is positive, then out_dir corresponds to a
    // ray exiting the surface into vaccum (ior = 1). However, note that
    // you should actually treat this case as _entering_ the surface, because
    // you want to compute the 'input' direction that would cause this output,
    // and to do so you can simply find the direction that out_dir would refract
    // _to_, as refraction is symmetric.
    Vec3 normal(0.0f, 1.0f, 0.0f);
    float cos_theta_i = dot(out_dir, normal);
    if (cos_theta_i < 0.0f)
        normal = -normal;

    float eta_i = cos_theta_i > 0 ? index_of_refraction : 1.0f;
    float eta_t = cos_theta_i > 0 ? 1.0f : index_of_refraction;
    float eta = eta_i / eta_t;

    cos_theta_i = abs(cos_theta_i);

    float sin_theta_t2 = eta * eta * (1 - cos_theta_i * cos_theta_i);
    if (sin_theta_t2 >= 1) {
        was_internal = false;
        return Vec3(0, 0, 0);
    }

    float cos_theta_t = sqrt(1 - sin_theta_t2);
    Vec3 refracted_dir = eta * -out_dir + (eta * cos_theta_i - cos_theta_t) * normal;
    
    was_internal = true;
    return refracted_dir;
}

BSDF_Sample BSDF_Lambertian::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 5
    // Implement lambertian BSDF. Use of BSDF_Lambertian::sampler may be useful

    BSDF_Sample ret;
    float pdf;
    ret.direction = sampler.sample(pdf);       // What direction should we sample incoming light from?
    ret.attenuation = evaluate(out_dir, ret.direction);        // What is the ratio of reflected/incoming light?
    ret.pdf = pdf;               // Was was the PDF of the sampled direction?
    return ret;
}

Spectrum BSDF_Lambertian::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    return albedo * (1.0f / PI_F);
}

BSDF_Sample BSDF_Mirror::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6
    // Implement mirror BSDF

    BSDF_Sample ret;
    ret.attenuation = reflectance; // What is the ratio of reflected/incoming light?
    ret.direction = reflect(out_dir);       // What direction should we sample incoming light from?
    ret.pdf = 1.0f; // Was was the PDF of the sampled direction? (In this case, the PMF)
    return ret;
}

Spectrum BSDF_Mirror::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // Technically, we would return the proper reflectance
    // if in_dir was the perfectly reflected out_dir, but given
    // that we assume these are single exact directions in a
    // continuous space, just assume that we never hit them
    // _exactly_ and always return 0.
    return {};
}

BSDF_Sample BSDF_Glass::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6
    BSDF_Sample ret;
    ret.attenuation = Spectrum(); // What is the ratio of reflected/incoming light?
    ret.direction = Vec3();       // What direction should we sample incoming light from?
    ret.pdf = 0.0f; // Was was the PDF of the sampled direction? (In this case, the PMF)

    // Implement glass BSDF.
    // (1) Compute Fresnel coefficient. Tip: use Schlick's approximation.
    auto schlick_fresnel_approx = [](float cos_theta_i, float n1, float n2) {
        if (cos_theta_i > 0.0)
            std::swap(n1, n2);

        float r0 = (n1 - n2) / (n1 + n2);
        r0 = r0 * r0;
        float x = 1.0f - cos_theta_i;
        return r0 + (1.0f - r0) * x * x * x * x * x;
    };

    Vec3 normal(0, 1, 0);
    float cos_theta_i = dot(out_dir, normal);
    float eta_i = 1.0f;
    float eta_t = index_of_refraction;
    if (cos_theta_i < 0) {
        // Ray is entering the medium
        cos_theta_i = -cos_theta_i;
        std::swap(eta_i, eta_t);
        normal = -normal;
    }

    float schlick_fresnel = schlick_fresnel_approx(cos_theta_i, index_of_refraction, 1.0f);
    
    // (2) Reflect or refract probabilistically based on Fresnel coefficient. Tip: RNG::coin_flip
    if (RNG::coin_flip(schlick_fresnel)) {
        // Reflect
        ret.direction = reflect(out_dir);
        ret.attenuation = reflectance * schlick_fresnel;
        ret.pdf = schlick_fresnel;
    } else {
        // Refract
        bool was_internal = false;
        ret.direction = refract(out_dir, index_of_refraction, was_internal);
        ret.attenuation = transmittance;
        if (was_internal) {
            // Total internal reflection occurred, reflect instead
            ret.direction = reflect(out_dir);
            ret.attenuation = Spectrum(1.0f, 1.0f, 1.0f);
            ret.pdf = schlick_fresnel;
        } else {
            ret.pdf = 1.0f - schlick_fresnel;
        }   
        
    }
    // (3) Compute attenuation based on reflectance or transmittance
    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
    return ret;
}

Spectrum BSDF_Glass::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // As with BSDF_Mirror, just assume that we never hit the correct
    // directions _exactly_ and always return 0.
    return {};
}

BSDF_Sample BSDF_Diffuse::sample(Vec3 out_dir) const {
    BSDF_Sample ret;
    ret.direction = sampler.sample(ret.pdf);
    ret.emissive = radiance;
    ret.attenuation = {};
    return ret;
}

Spectrum BSDF_Diffuse::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // No incoming light is reflected; only emitted
    return {};
}

BSDF_Sample BSDF_Refract::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6
    // Implement pure refraction BSDF.

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?

    BSDF_Sample ret;
    bool was_internal = false;

    // Figure out which $\eta$ is incident and which is transmitted

    ret.direction = refract(out_dir, index_of_refraction, was_internal); 
    
    // Compute ray direction for specular transmission
    if (was_internal){
        ret.attenuation = Spectrum(1.0f, 1.0f, 1.0f);
        ret.direction = reflect(out_dir);
    } else {
        ret.attenuation = transmittance;
        //ret.attenuation = ret.attenuation / abs(ret.direction.y);
    }
    
    return ret;
}

Spectrum BSDF_Refract::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // As with BSDF_Mirror, just assume that we never hit the correct
    // directions _exactly_ and always return 0.
    return {};
}

} // namespace PT
