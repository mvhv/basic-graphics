// Basic C++ Implementation of Backwards Ray-tracing
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <sstream>

#define MAX_RAY_DEPTH 5


// blends doubles with the specified ratio
double blendDoubles(double a, double b, double ratio) {
    return b * ratio + a * (1 - ratio);
}

// 3-vector representation for points and vectors in 3D space, as well as colours.
class Vec3 {
    public:
        double x, y, z;

        Vec3() {
            x = 0.;
            y = 0.;
            z = 0.;
        }

        Vec3(double xx, double yy, double zz) {
            x = xx;
            y = yy;
            z = zz;
        }

        Vec3(double weight) {
            x = weight;
            y = weight;
            z = weight;
        }

        Vec3 operator * (const double multiplier) const {
            return Vec3(x * multiplier, y * multiplier, z * multiplier);
        }

        Vec3 operator * (const Vec3 &vec) const {
            return Vec3(x * vec.x, y * vec.y, z * vec.z);
        }

        Vec3 operator - (const Vec3 &vec) const {
            return Vec3(x - vec.x, y - vec.y, z - vec.z);
        }

        Vec3 operator + (const Vec3 &vec) const {
            return Vec3(x + vec.x, y + vec.y, z + vec.z);
        }

        Vec3& operator += (const Vec3 &vec) {
            x += vec.x;
            y += vec.y;
            z += vec.z;
            return *this;
        }

        Vec3& operator *= (const Vec3 &vec) {
            x *= vec.x;
            y *= vec.y;
            z *= vec.z;
            return *this;
        }

        // urnary negative
        Vec3 operator - () const {
            return Vec3(-x, -y, -z);
        }

        // calculate the dot product between this and vec
        double dot (const Vec3 &vec) const {
            return (x * vec.x + y * vec.y + z * vec.z);
        }

        // calulate the square of the length
        double length2() const {
            return (x * x + y * y + z * z);
        }

        double length() const {
            return std::sqrt(length2());
        }

    // DO SOMETHING ------------------------------------------------------------------------------------
        friend std::ostream & operator << (std::ostream &os, const Vec3 &v) {
            os << "[" << v.x << " " << v.y << " " << v.z << "]";
            return os;
        }

        // normalise vector fields so that length == 1
        Vec3& normalise() {
            if (length2() > 0) {
                double normalisationFactor = 1 / length();
                x *= normalisationFactor;
                y *= normalisationFactor;
                z *= normalisationFactor;
            }
            return *this;
        }
};


class Ray {
    public:
        Vec3 origin, direction;

        Ray(const Vec3 direction_, const Vec3 origin_ = Vec3(0., 0., 0.)) {
            origin = origin_;
            direction = direction_;
        }
};


class Sphere {
    public:
        Vec3 centre, surfaceColour, emissionColour; // centre position, surface colour and emissions
        double radius, radius2, transparency, reflectivity; // radius and radius squared, also surface transparency and reflectivity

        Sphere(const Vec3& centre_, const double& radius_, const Vec3& surfaceColour_, const double& reflectivity_,
                const double& transparency_, const Vec3& emissionColour_) {
            centre = centre_;
            radius = radius_;
            radius2 = radius_ * radius_;
            surfaceColour = surfaceColour_;
            emissionColour = emissionColour_;
            transparency = transparemcy_;
            reflectivity = reflectivity_;
        }

        // compute ray-sphere intersection using simple geometric solution
        centreBearing = centre - rayOrigin
        timeCentreAdjacent = centerBearing . rayDirection
        if tca < 0 then sphere is behind the origin

        
        bool intersection(const Ray& ray, float &t0, float &t1) const {
            
            Vec3 centreBearing = centre - Ray.origin;
            double Tca, 




            Vec3f l = centre - rayorig;
            float tca = l.dot(raydir);
            if (tca < 0) return false;
            float d2 = l.dot(l) - tca * tca;
            if (d2 > radius2) return false;
            float thc = std::sqrt(radius2 - d2);
            t0 = tca - thc;
            t1 = tca + thc;
            return true;
        }
};


Vec3f trace(
    const Vec3f &rayorig,
    const Vec3f &raydir,
    const std::vector<Sphere> &spheres,
    const int &depth)
{
    // if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    float tnear = INFINITY;
    const Sphere* sphere = nullptr;
    // find intersection of this ray with the sphere in the scene
    for (unsigned i = 0; i < spheres.size(); ++i) {
        float t0 = INFINITY, t1 = INFINITY;
        if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &spheres[i];
            }
        }
    }

    //if there's no intersection return black or background colour
    if (!sphere) return Vec3f(2);
    Vec3f surfaceColour = 0; // colour of the ray/surface of the object intersected by the ray
    Vec3f phit = rayorig + raydir * tnear; // point of intersection
    Vec3f nhit = phit - sphere->centre; // normal at the intersection point
    nhit.normalise(); // normalise normal direction

    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want positive.
    float bias = 1e-4; //add some bias to the point from which we will be tracing
    bool inside = false;

    if (raydir.dot(nhit) > 0) {
    nhit = -nhit;
    inside = true;
    }
    if ((sphere->transparency > 0 || sphere->reflectivity > 0) && depth < MAX_RAY_DEPTH) {
        float facingRatio = -raydir.dot(nhit);
        // change the mix value to tweak the effect
        float fresnelEffect = blendDoubles(pow(1 - facingRatio, 3), 1, 0.1);
        // compute reflection direction (not need to normalise because all vectors
        // are already normalised)
        Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit); 
        refldir.normalise();
        Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
        Vec3f refraction = 0;
        // if the sphere is also transparent compute refraction ray (transmission)
        if (sphere->transparency) {
            float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
            float cosi = -nhit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vec3f refrdir = raydir * eta + nhit * (eta * cosi - std::sqrt(k));
            refrdir.normalise();
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
        }

        // the result is a mix of reflection and refraction (if the sphere is transparent)
        surfaceColour = (reflection * fresnelEffect + refraction * (1 - fresnelEffect) * sphere->transparency) * sphere->surfaceColour;
    } else {
        // it's a diffuse object, no need to raytrace any further
        for (unsigned i = 0; i < spheres.size(); ++i) {
            if (spheres[i].emissionColour.x > 0) {
                // this is a light
                Vec3f transmission = 1;
                Vec3f lightDirection = spheres[i].centre - phit;
                lightDirection.normalise();
                for (unsigned j = 0; j < spheres.size(); ++j) {
                    if (i != j) {
                        float t0, t1;
                        if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                            transmission = 0;
                            break;
                        }
                    }
                }
                surfaceColour += sphere->surfaceColour * transmission * std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColour;
            }
        }
    }

    return surfaceColour + sphere->emissionColour;
}

long floatToStringClamped(const float num) {
    return std::lrint(std::min(255.f, num * 255.f));
}


std::string ppmFormat(const Vec3f &image) {
    std::stringstream ss;
    ss << floatToStringClamped(image.x) << " " << floatToStringClamped(image.y) << " " << floatToStringClamped(image.z) << "\n";
    return ss.str();
}

// main rendering function we compute a camera ray for each pixel of the image, trace it and return a colour.
// If the ray hits a sphere, we return the colour of the sphere at the intersection point, else we return the background colour.

void render(const std::vector<Sphere> &spheres) {
    unsigned width = 640, height = 480;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectRatio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // trace rays

    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            // pixel center points
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectRatio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            raydir.normalise();
            *pixel = trace(Vec3f(0), raydir, spheres, 0);
        }
    }

    // Save result to a PPM image ( keep these flags if you compiule under Windows)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    ofs << "P3\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << ppmFormat(image[i]);
    }
    ofs.close();
    delete [] image;
}

// in the main function, we will create the scene which is composed of 5 spheres and 1 light (which is also a sphere).
// Then, once the scene description is complete we render that scene, by calling the render() function.
int main(int argc, char **argv) {
    std::vector<Sphere> spheres;
    // position, radius, surface, colour, reflectivity, transparency, emission colour
    spheres.push_back(Sphere(Vec3f(0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0)); 
    spheres.push_back(Sphere(Vec3f(0.0, 0, -20), 4, Vec3f(1.00, 0.32, 0.36), 1, 0.5)); 
    spheres.push_back(Sphere(Vec3f(5.0, -1, -15), 2, Vec3f(0.90, 0.76, 0.46), 1, 0.0)); 
    spheres.push_back(Sphere(Vec3f(5.0, 0, -25), 3, Vec3f(0.65, 0.77, 0.97), 1, 0.0)); 
    spheres.push_back(Sphere(Vec3f(-5.5, 0, -15), 3, Vec3f(0.90, 0.90, 0.90), 1, 0.0)); 
    // light
    spheres.push_back(Sphere(Vec3f(0.0, 20, -30), 3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
    render(spheres);
    return 0;
}