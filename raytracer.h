#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

#define screen_width 128
#define screen_height 128

// Vector3D same as struct
struct Vector3D
{
    double x, y, z;
};

struct RGB
{
    double r, g, b;
};

struct Triangle
{
    Vector3D A, B, C;
};

struct Radiance3
{
    double red, green, blue;
};


struct Light
{
    Vector3D position;
    Radiance3 colour;
};


struct Material
{
    Radiance3 emission;
    Radiance3 diffuse;
    Radiance3 specular;
};

struct AreaLight
{
    Triangle triangle;
    Radiance3 colour;
    Material lightMaterial;
};

struct Surface
{
    Triangle triangle;
    Material material;
    double extinctionProbability;
    bool reflects;
    bool emits;
};

struct Surfel
{
    Vector3D position;
    Material material;
    Vector3D normal;
    double extinctionProbability;
    bool reflects;
    bool emits;
    bool valid;

};

struct Ray
{
    Vector3D rayDirection;
    Vector3D origin;
};


Radiance3 screen[screen_width][screen_height];
Vector3D camera;
Vector3D light;
std::vector<Surface> surfaceVector;
std::vector<Light> pointLightVector;
std::vector<AreaLight> areaLightVector;

Radiance3 pathTrace(Ray ray, bool isEyeRay);
void createWorld();
void initialise();
Radiance3 clamp(Radiance3 temp);
void exportFile(std::string outname);
Vector3D convertPCS(Vector3D p, Vector3D u, Vector3D w, Vector3D n, Vector3D o);
bool Barycentric(Vector3D oPCS, Vector3D pNCS, Vector3D qNCS, Vector3D rNCS);
Surfel intersect(Ray ray);
bool lineOfSight(Surfel surfel, Surfel lightSurfel);
Radiance3 estimateDirectPointLight(Surfel surfel, Ray ray);
Radiance3 estimateDirectAreaLight(Surfel surfel, Ray ray);
Radiance3 BSDF(Surfel surfel, Vector3D omega, Vector3D rayDirection);
Radiance3 estimateIndirectLight(Surfel surfel, Ray ray, bool isEyeRay);

Vector3D operator - (Vector3D c1, Vector3D c2)
{
    Vector3D out = {c1.x - c2.x, c1.y - c2.y, c1.z - c2.z};
    return out;
}

Vector3D operator + (Vector3D c1, Vector3D c2)
{
    Vector3D out = {c1.x + c2.x, c1.y + c2.y, c1.z + c2.z};
    return out;
}

Vector3D operator / (Vector3D vector, double divisor)
{
    Vector3D out = {vector.x / divisor, vector.y / divisor, vector.z / divisor};
    return out;
}

Vector3D operator * (Vector3D vector, double doub)
{
    Vector3D out = {vector.x * doub, vector.y * doub, vector.z * doub};
    return out;
}

Vector3D operator * (double doub, Vector3D vector)
{
    Vector3D out = {vector.x * doub, vector.y * doub, vector.z * doub};
    return out;
}



Radiance3 operator += (Radiance3 r1, Radiance3 r2)
{
    Radiance3 out;
    double r = r1.red + r2.red;
    double g = r1.green + r2.green;
    double b = r1.blue + r2.blue;
    out = {r, g, b};
    return out;
}

Radiance3 operator + (Radiance3 r1, Radiance3 r2)
{
    Radiance3 out = {r1.red + r2.red, r1.green + r2.green, r1.blue + r2.blue};
    return out;
}

// Function to simplify dot product
double DotProduct(Vector3D pointA, Vector3D pointB)
{
    double result = 0;
    result += pointA.x * pointB.x;
    result += pointA.y * pointB.y;
    result += pointA.z * pointB.z;
    return result;
}

double VectorLength(Vector3D point)
{
    double result = sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
    return result;
}

double max(double num1, double num2)
{

    if(num1 > num2)
    {
        return num1;
    }
    else if (num2 > num1)
    {
        return num2;
    }
    else
        return num1;
}

double randomDouble()
{

    return rand() / (RAND_MAX + 1.0);
}



// Function to simplify cross product
Vector3D CrossProduct(Vector3D vectorA, Vector3D vectorB)
{
    Vector3D result = {0, 0, 0};
    result.x = (vectorA.y * vectorB.z) - (vectorA.z * vectorB.y);
    result.y = (vectorA.z * vectorB.x) - (vectorA.x * vectorB.z);
    result.z = (vectorA.x * vectorB.y) - (vectorA.y * vectorB.x);
    return result;
}

Vector3D normalise(Vector3D vector)
{
    double length = sqrt((vector.x*vector.x) + (vector.y*vector.y) + (vector.z * vector.z));
    vector.x = vector.x/length;
    vector.y = vector.y/length;
    vector.z = vector.z/length;
    return vector;
}

void vectorPrint(Vector3D vector)
{
    std::cout << "vector x:" << vector.x << " y:" << vector.y << " z:" << vector.z << '\n';
}

void radiancePrint(Radiance3 rad)
{
    std::cout << "radiance red:" << rad.red << " green:" << rad.green << " blue:" << rad.blue << '\n';
}


Vector3D invert(Vector3D vector)
{
    vector.x = - vector.x;
    vector.y = - vector.y;
    vector.z = - vector.z;
    return vector;
}

Radiance3 operator / (Radiance3 rad3, double divisor)
{
    Radiance3 output;
    output.red = rad3.red/divisor;
    output.green = rad3.green/divisor;
    output.blue = rad3.blue/divisor;
    return output;
}

Radiance3 operator * (Radiance3 rad3, double multiplier)
{
    Radiance3 output;
    output.red = rad3.red * multiplier;
    output.green = rad3.green * multiplier;
    output.blue = rad3.blue * multiplier;
    return output;
}

Radiance3 operator * (Radiance3 rad1, Radiance3 rad2)
{
    Radiance3 output;
    output.red = rad1.red * rad2.red;
    output.green = rad1.green * rad2.green;
    output.blue = rad1.blue * rad2.blue;
    return output;
}
