#include "raytracer.h"
Material lightMaterial;
#define nSamples 100


int main(int argc, char **argv)
{
    // Position camera/eye at origin
    camera = {0, 0, 0};

    // Seed rand function for use later with area lights
    srand((unsigned)time(NULL));
    // Call function to build world by defining triangles
    createWorld();

    // Initialise screen values all to 0
    initialise();

    // Loop through every pixel on the screen
    for(int y = 0; y < screen_height; y++)
    {
        for(int x = 0; x < screen_width; x++)
        {
            // Initialise origin to 0
            Vector3D origin = {0, 0, 0};
            // Offset so that field of view is looking into the box with all sides visible
            Vector3D rayVector = { ( ((((double)screen_width-x)/(double)screen_width) - 0.5) - camera.x), (((((double)screen_width-y)/(double)screen_width) - 0.5) - camera.y), (0.5) };
            //std::cout << "ray x: " << rayVector.x << " y: " << rayVector.y << " z: " << rayVector.z << '\n';
            // Create Ray from origin
            Ray ray = {rayVector, origin};
            Radiance3 pixelRadiance = {0, 0, 0};

            // For the specified number of samples, call the pathTrace function to calculate lighting
            for(int sample = 0; sample < nSamples; sample++)
            {
                pixelRadiance = pixelRadiance + pathTrace(ray, true);
                //std::cout << "R: " << pixelRadiance.red << " G: " << pixelRadiance.green << " B: " << pixelRadiance.blue << '\n';

                // Set pixel value to calculated RGB
                screen[x][y] = pixelRadiance;
            }
        }
    }
    // Output the screen array to a ppm file
    exportFile("box.ppm");
}

// Function to initialise all screen values to 0
void initialise()
{
    Radiance3 Default = {0, 0, 0};
    for(int y = 0; y < screen_height; y++)
    {
        for(int x = 0; x < screen_width; x++)
        {
            screen[x][y] = Default;
        }

    }
}

Radiance3 pathTrace(Ray ray, bool isEyeRay)
{
    Surfel surfel;

    // Initial radiance output of 0
    Radiance3 output = {0, 0, 0};
    //std::cout << "pre intersect" << '\n';

    // Store what the ray intersects in surfel
    surfel = intersect(ray);
    if(surfel.valid)
    {

        // If ray hits light, add emission
        if(isEyeRay && surfel.emits)
        {
            //std::cout << "emits" << '\n';
            output = output + surfel.material.emission;
        }

        // Calculate direct lighting from point and area lights
        if((!isEyeRay) || surfel.reflects)
        {
            output = output + estimateDirectPointLight(surfel, ray);
            output = output + estimateDirectAreaLight(surfel, ray);
            //radiancePrint(output);

        }

        // Add indirect illumination
        if(!(isEyeRay) || true)
        {
            output = output + estimateIndirectLight(surfel, ray, isEyeRay);
        }

    }

    //std::cout << "after intersect" << '\n';

    // Return final radiance for pixel
    return output;

}

/*  This functional calculates the first surfel a ray intersects by looping through
*   all surfaces and lights and testing for intersections. It then determines which
*   surfel intersected is closest and returns it
*/
Surfel intersect(Ray ray)
{
    std::vector<Surfel> intersectVector;
    std::vector<double> distanceVector;
    Surfel surfel;

    // Test for intersection with any surfaces
    for(int i=0; i < surfaceVector.size(); i++)
    {
        // Construct vectors U, V and W along with the normal making sure to normalise
        Triangle triangle = surfaceVector[i].triangle;
        Vector3D u = normalise(triangle.B - triangle.A);
        Vector3D v = normalise(triangle.C - triangle.A);
        Vector3D normal = normalise(CrossProduct(u, v));
        Vector3D w = normalise(CrossProduct(u, normal));

        //Find Point O before converting to PCS
        //s + l * lscale is formula used to find pointO
        Vector3D pointO = ray.origin + ray.rayDirection*(DotProduct(triangle.A, normal)/DotProduct(ray.rayDirection, normal));
        Vector3D oPCS = convertPCS(triangle.A, u, w, normal, pointO);
        Vector3D pPCS = convertPCS(triangle.A, u, w, normal, triangle.A);
        Vector3D qPCS = convertPCS(triangle.A, u, w, normal, triangle.B);
        Vector3D rPCS = convertPCS(triangle.A, u, w, normal, triangle.C);

        //std::cout << "point O x: " << pointO.x << " y: " << pointO.y << " z: " << pointO.z << '\n';

        // Use Barycentric to test whether an intersection
        if(Barycentric(oPCS, pPCS, qPCS, rPCS))
        {
            // If intersection then make sure no self intersection before storing surfel

            //std::cout << "intersect" << '\n';
            surfel = {pointO, surfaceVector[i].material, invert(normal), surfaceVector[i].extinctionProbability, surfaceVector[i].reflects, surfaceVector[i].emits, true};
            double length =  VectorLength(pointO - ray.origin);
            // make sure no self intersection
            if(length > 0.1)
            {
                intersectVector.push_back(surfel);
                distanceVector.push_back(length);
            }

        }
    }
    // Test for intersection with any lights
    for(int i=0; i < areaLightVector.size(); i++)
    {
        // Construct vectors U, V and W along with the normal making sure to normalise
        Triangle triangle = areaLightVector[i].triangle;
        Vector3D u = normalise(triangle.B - triangle.A);
        Vector3D v = normalise(triangle.C - triangle.A);
        Vector3D normal = normalise(CrossProduct(u, v));
        Vector3D w = normalise(CrossProduct(u, normal));

        //Find Point O before converting to PCS
        //s + l * lscale
        Vector3D pointO = ray.origin + ray.rayDirection*(DotProduct(triangle.A, normal)/DotProduct(ray.rayDirection, normal));
        Vector3D oPCS = convertPCS(triangle.A, u, w, normal, pointO);
        Vector3D pPCS = convertPCS(triangle.A, u, w, normal, triangle.A);
        Vector3D qPCS = convertPCS(triangle.A, u, w, normal, triangle.B);
        Vector3D rPCS = convertPCS(triangle.A, u, w, normal, triangle.C);

        //std::cout << "point O x: " << pointO.x << " y: " << pointO.y << " z: " << pointO.z << '\n';

        // Use Barycentric to test for intersection
        if(Barycentric(oPCS, pPCS, qPCS, rPCS))
        {
            // If intersection then make sure no self intersection before storing surfel

            //std::cout << "light intersect" << '\n';
            surfel = {pointO, areaLightVector[i].lightMaterial, invert(normal), 50, true, true, true};
            double length =  VectorLength(pointO - ray.origin);
            // make sure no self intersection
            if(length > 0)
            {
                intersectVector.push_back(surfel);
                distanceVector.push_back(length);
            }

        }
    }
    // If empty then no intersection so return empty surfel
    if(distanceVector.empty())
    {
        Surfel blank;
        blank.valid = false;
        return blank;
        //throw "no intersection";
    }

    // Return closest intersected surfel
    double dist = 99999999;
    int index;
    for(int j = 0; j < distanceVector.size(); j++)
    {
        if(distanceVector[j] < dist)
        {
            dist = distanceVector[j];
            index = j;
        }
    }
    return intersectVector[index];
}

bool lineOfSight(Surfel surfel, Surfel lightSurfel)
{
    Ray ray = {lightSurfel.position - surfel.position, surfel.position};
    ray.rayDirection = normalise(ray.rayDirection);
    //std::cout << "ray x: " << ray.rayDirection.x << " y: " << ray.rayDirection.y << " z: " << ray.rayDirection.z << '\n';
    double length =  VectorLength(ray.rayDirection);

    for(int i=0; i < surfaceVector.size(); i++)
    {
        // Construct vectors U, V and W along with the normal making sure to normalise
        Triangle triangle = surfaceVector[i].triangle;
        Vector3D u = normalise(triangle.B - triangle.A);
        Vector3D v = normalise(triangle.C - triangle.A);
        Vector3D normal = normalise(CrossProduct(u, v));
        Vector3D w = normalise(CrossProduct(u, normal));

        //Find Point O before converting all to PCS
        //s + l * lscale
        Vector3D pointO = ray.origin + ray.rayDirection*(DotProduct(triangle.A, normal)/DotProduct(ray.rayDirection, normal));
        Vector3D oPCS = convertPCS(triangle.A, u, w, normal, pointO);
        Vector3D pPCS = convertPCS(triangle.A, u, w, normal, triangle.A);
        Vector3D qPCS = convertPCS(triangle.A, u, w, normal, triangle.B);
        Vector3D rPCS = convertPCS(triangle.A, u, w, normal, triangle.C);

        // Use Barycentric to test for intersection
        if(Barycentric(oPCS, pPCS, qPCS, rPCS) )
        {
            // Test whether intersected object is between the lightsurfel and surfel

            //Distance between surface and object
            double Light2Object =  VectorLength(lightSurfel.position - pointO );
            //distance between surfel and object
            double Surfel2Object = VectorLength(surfel.position - pointO);
            if(Surfel2Object > 1 && Light2Object > 1)
            {
                //std::cout << "Light2Object " << Light2Object << '\n';
                //std::cout << "Surfel2Object " << Surfel2Object << '\n';
                //std::cout << "length " << length << '\n';

                // If blocked return false so shadow
                if(Light2Object < (length) && Surfel2Object < (length))
                {
                    //std::cout << "false\n" << '\n';
                    return false;
                }
            }
        }
    }
    //std::cout << "true\n" << '\n';
    return true;
}

// Calculates the lighting from point light sources for a pixel
Radiance3 estimateDirectPointLight(Surfel surfel, Ray ray)
{
    Radiance3 output = {0, 0, 0};
    Surfel check;

    // Loop through all point lights for lighting
    for(int source = 0; source < pointLightVector.size(); source++)
    {
        //std::cout << "light" << '\n';
        // Cast shadow ray to test visibility
        Ray shadowRay = {pointLightVector[source].position - surfel.position, surfel.position};

        check = intersect(shadowRay);
        if(surfel.valid)
        {
            // Calculate incoming vector and BSDF to convert to light output if not in shadow

            Vector3D omega_i = pointLightVector[source].position - surfel.position;
            double distance = VectorLength(omega_i);
            Radiance3 E_i = (pointLightVector[source].colour / (4 * M_PI * distance * distance )) *500000;
            //std::cout << "source r: " <<E_i.red << " g: " << E_i.green << " b: " << E_i.blue << '\n';
            //Radiance3 bsdftemp = BSDF(surfel, omega_i, invert(ray.rayDirection));
            //std::cout << "bsdf r: " << bsdftemp.red << " g: " << bsdftemp.green << " b: " << bsdftemp.blue << '\n';
            output = output + BSDF(surfel, omega_i, invert(ray.rayDirection)) * E_i * max(0, DotProduct(normalise(omega_i), surfel.normal));

        }
    }

    return output;
}


// Calculate the lighting from area lights
Radiance3 estimateDirectAreaLight(Surfel surfel, Ray ray)
{
    Radiance3 output = {0, 0, 0};
    Surfel check;

    // Loop through all area lights
    for(int source = 0; source < areaLightVector.size(); source++)
    {
        //Ray shadowRay = {areaLightVector[source].position - surfel.position, surfel.position};
        // Choose random point on light
        Vector3D pointA = areaLightVector[source].triangle.A;
        Vector3D pointB = areaLightVector[source].triangle.B;
        Vector3D pointC = areaLightVector[source].triangle.C;

        // x = v1 + a(v2 - v1) + b(v3 - v1) to choose random point
        Vector3D pointRandom = pointA + (randomDouble() * (pointB - pointA)) + (randomDouble() * (pointC - pointA));

        // Calulate vectors u, v and normal for the random surfel on light chosen
        Vector3D u = normalise(pointB - pointA);
        Vector3D v = normalise(pointC - pointA);
        Vector3D normal = {0, -1, 0};
        normal = normalise(CrossProduct(u, v));
        Surfel lightSurfel = {pointRandom, lightMaterial, invert(normal), 50, false, true, true };


        // Cast shadow ray and test for intersection
        Ray shadowRay = {lightSurfel.position - surfel.position, surfel.position};
        check = intersect(shadowRay);

        //std::cout << "light surfel x: " << lightSurfel.position.x << " y: " << lightSurfel.position.y << " z: " << lightSurfel.position.z << '\n';
        //std::cout << "surfel x: " << surfel.position.x << " y: " << surfel.position.y << " z: " << surfel.position.z << '\n';

        // If in line of sight without shadow, calculate the incoming vector and BSDF to determine radiance
        if(lineOfSight(surfel, lightSurfel))
        {
            Vector3D omega_i = lightSurfel.position - surfel.position;
            //std::cout << "omega_i x: " << omega_i.x << " y: " << omega_i.y << " z: " << omega_i.z << '\n';
            double distance = VectorLength(omega_i);

            Radiance3 bsdftemp;
            bsdftemp = BSDF(surfel, omega_i, invert(ray.rayDirection));

            //std::cout << "bsdf r: " << bsdftemp.red << " g: " << bsdftemp.green << " b: " << bsdftemp.blue << '\n';
            //std::cout << "light dot product " << DotProduct(normalise(omega_i), lightSurfel.normal) << '\n';

            output = output + clamp(bsdftemp) * 0.1 * M_PI * max(0, DotProduct(normalise(omega_i), lightSurfel.normal));
        }

    }
    //std::cout << "output:" << '\n';
    //radiancePrint(output);
    //std::cout << " " << '\n';
    return output;

}

// Function for clamping any negative values to 0 to prevent any errors stopping output image opening
Radiance3 clamp(Radiance3 temp)
{
    if(temp.red < 0)
    {
        temp.red = 0;
    }
    if(temp.green < 0)
    {
        temp.green = 0;
    }
    if(temp.blue < 0)
    {
        temp.blue = 0;
    }

    return temp;
}

// Calculate any indirect lighting for pixel
Radiance3 estimateIndirectLight(Surfel surfel, Ray ray, bool isEyeRay)
{
    Radiance3 output = {0, 0, 0};
    Vector3D bounceVector;

    // Generate random number, if higher than extinction probability then extinguish photon
    if(rand() % 100 > surfel.extinctionProbability)
    {
        //std::cout << "extinct" << '\n';
        return output;
    }
    else
    {   // Not extinguished
        /* Generate random vector and use dot product to test whether in valid direction
        *  and if not, generate random vectors until one is
        */
        while(true)
        {
            bounceVector = {(rand() % 100) - 50.0, (rand() % 100) - 50.0, (rand() % 100) - 50.0};
            bounceVector = normalise(bounceVector);
            //std::cout << "bounce vector x: " << bounceVector.x << " y:" << bounceVector.y << " z: " << bounceVector.z <<'\n';
            //std::cout << "normal" << '\n';
            //vectorPrint(surfel.normal);
            //std::cout << "dot product: " << DotProduct(surfel.normal, bounceVector) << '\n';
            if(DotProduct(surfel.normal, bounceVector) > 0)
            {
                break;
            }
        }
        //std::cout << "surfel position" << '\n';
        //vectorPrint(surfel.position);

        // Return the ray that has bounced to pathTrace but with a false as its no longer the eye ray
        Ray bounceRay = {bounceVector, surfel.position};
        return pathTrace(bounceRay, false);
    }
}

Radiance3 BSDF(Surfel surfel, Vector3D omega, Vector3D rayDirection)
{
    //diffuse + specular for phong
    //omega is Vl, rayDirection is Ve

    // Calculate diffuse lighting
    double nv = DotProduct(surfel.normal, omega);
    double normalLength = VectorLength(surfel.normal);
    double vlLength = VectorLength(omega);

    //Radiance3 =  (light.colour)lightBrightness * (surfel.material.diffuse)materialBrightness * (nv/normalLength * vlLength)
    Radiance3 diffuse = surfel.material.diffuse * (nv/(normalLength * vlLength));
    //std::cout << "diffuse r: " << diffuse.red << " g: " << diffuse.green << " b: " << diffuse.blue << '\n';

    // Calculate specular lighting
    // Vl = omega
    // Ve = e - p
    Vector3D Vb = (omega + rayDirection)/2.0;
    double nVb = DotProduct(surfel.normal, Vb);
    double VbLength = VectorLength(Vb);
    Radiance3 specular = surfel.material.specular * pow((nVb / (normalLength * VbLength)), 20);
    //std::cout << "specular r: " << specular.red << " g: " << specular.green << " b: " << specular.blue << '\n';
    return (diffuse+specular);


}


// Function that uses barycentric coordinates to test whether point O intersects triangle
bool Barycentric(Vector3D oPCS, Vector3D pNCS, Vector3D qNCS, Vector3D rNCS)
{
    // B - A
    double gradh0 = qNCS.x - pNCS.x;
    double grad0 = - (qNCS.y - pNCS.y);
    double c0 = (grad0*pNCS.x) + (gradh0*pNCS.y);

    // C - B
    double gradh1 = rNCS.x - qNCS.x;
    double grad1 = -(rNCS.y - qNCS.y);
    double c1 = (grad1*qNCS.x) + (gradh1*qNCS.y);

    // A - C
    double gradh2 = pNCS.x - rNCS.x;
    double grad2 = -(pNCS.y - rNCS.y);
    double c2 = (grad2*rNCS.x) + (gradh2*rNCS.y);

    // Calculate the distance from each point to the line opposite in the triangel
    double CtoBA = (( ((grad0 * rNCS.x) + (gradh0 * rNCS.y)) - c0 ) / (sqrt(grad0*grad0 + gradh0*gradh0)) );
    double AtoCB = (((grad1 * pNCS.x) + (gradh1 * pNCS.y) - c1 ) / sqrt(grad1*grad1 + gradh1*gradh1));
    double BtoAC = (((grad2 * qNCS.x) + (gradh2 * qNCS.y) - c2 ) / sqrt(grad2*grad2 + gradh2*gradh2));

    // Calculate distance from line to P
    double BA = (( ((grad0 * oPCS.x) + (gradh0 * oPCS.y)) - c0 ) / (sqrt(grad0*grad0 + gradh0*gradh0)) );
    double CB = (((grad1 * oPCS.x) + (gradh1 * oPCS.y) - c1 ) / sqrt(grad1*grad1 + gradh1*gradh1));
    double AC = (((grad2 * oPCS.x) + (gradh2 * oPCS.y) - c2 ) / sqrt(grad2*grad2 + gradh2*gradh2));

    // Calculate alpha, beta, gamma using distances
    double alpha = CB / AtoCB;
    double beta =  AC / BtoAC;
    double gamma = BA / CtoBA;

    // If true then intersects, else false
    if( (alpha > 0) && (beta > 0) && (gamma > 0) )
    {
        return true;
    }
    else
    return false;

}

// FUnction that converts coordinates onto planar coordinate system using dot product
Vector3D convertPCS(Vector3D p, Vector3D u, Vector3D w, Vector3D n, Vector3D o)
{
    Vector3D s = o-p;
    Vector3D s_PCS = {DotProduct(s, u), DotProduct(s, w), DotProduct(s, n)};
    return s_PCS;
}

/* Large function used to define all triangles, surfaces and lights required to
*  build the scene.
*/
void createWorld()
{
    // create light
    Vector3D light1pos = {0, 0, 50};
    Radiance3 white = {1, 1, 1};
    Light light = {light1pos, white};
    //pointLightVector.push_back(light);

    Radiance3 noDiffuse = {0, 0, 0};
    Radiance3 noSpecular = {0, 0, 0};
    lightMaterial = {white, noDiffuse, noSpecular};

    // create triangle 1
    Vector3D A, B, C;
    Triangle triangle;
    Radiance3 emission, diffuse, specular;
    Material mat;
    Surface surface;
    //surfaceVector.push_back(surface);

    //create area lights
    A = {-15, 48.9, 70};
    B = {-15, 48.9, 50};
    C = {15, 48.9, 50};
    triangle = {A, B, C};
    AreaLight arealight = {triangle, white, lightMaterial};
    areaLightVector.push_back(arealight);


    A = {15, 48.9, 70};
    B = {-15, 48.9, 70};
    C = {15, 48.9, 50};
    triangle = {A, B, C};
    arealight = {triangle, white,lightMaterial};
    areaLightVector.push_back(arealight);

    // ceiling
    C = {-51, 49, 100};
    B = {-51, 49, 0};
    A = {51, 49, 0};
    triangle = {B, C, A};
    emission = {0, 0, 0};
    diffuse = {1, 1, 1};
    specular = {1, 1, 1};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    C = {51, 49, 100};
    B = {-51, 49, 100};
    A = {51, 49, 0};
    triangle = {A, B, C};
    emission = {0, 0, 0};
    diffuse = {1, 1, 1};
    specular = {1, 1, 1};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    // create ground triangles
    A = {-50, -50, 100};
    B = {-50, -50, 0};
    C = {50, -50, 0};
    triangle = {A, B, C};
    emission = {0, 0, 0};
    diffuse = {1, 1, 1};
    specular = {1, 1, 1};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    A = {50, -50, 100};
    B = {-51, -50, 100};
    C = {50, -50, 0};
    triangle = {A, B, C};
    emission = {0, 0, 0};
    diffuse = {1, 1, 1};
    specular = {1, 1, 1};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    //Create green right wall

    A = {-50, 50, 0};
    B = {-50, -51, 0};
    C = {-50, -51, 100};
    triangle = {A, B, C};
    emission = {0, 0, 0};
    diffuse = {0, 1, 0};
    specular = {0, 1, 0};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    A = {-50, 50, 0};
    B = {-50, -51, 100};
    C = {-50, 50, 100};
    triangle = {A, B, C};
    emission = {0, 0, 0};
    diffuse = {0, 1, 0};
    specular = {0, 1, 0};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    //Create red left wall

    C = {50, 50, 0};
    B = {50, -51, 0};
    A = {50, -51, 100};
    triangle = {A, B, C};
    emission = {0, 0, 0};
    diffuse = {1, 0, 0};
    specular = {1, 0, 0};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    C = {50, 50, 0};
    B = {50, -51, 100};
    A = {50, 50, 100};
    triangle = {A, B, C};
    emission = {0, 0, 0};
    diffuse = {1, 0, 0};
    specular = {1, 0, 0};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    // create white back wall
    A = {-51, 50, 99};
    B = {-51, -50, 99};
    C = {51, -51, 99};
    triangle = {A, B, C};
    emission = {0, 0, 0};
    diffuse = {1, 1, 1};
    specular = {1, 1, 1};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);

    A = {-51, 50, 99};
    B = {51, -51, 99};
    C = {51, 50, 99};
    triangle = {A, B, C};
    emission = {1, 1, 1};
    diffuse = {1, 1, 1};
    specular = {1, 1, 1};
    mat = {emission, diffuse, specular};
    surface = {triangle, mat, 50, true, false};
    surfaceVector.push_back(surface);


}


// Exports the screen as a PPM file whilst adding gamma correction
void exportFile(std::string outname)
{
    // Max and Min values used for adding gamma correction
    double redMax = 0;
    double greenMax = 0;
    double blueMax = 0;

    double redMin = 999;
    double greenMin = 999;
    double blueMin = 999;

    for(int j = 0; j<screen_width; j++)
    {
        // Print each value on the row
        for(int i = 0; i<screen_height; i++)
        {
            if(screen[i][j].red > redMax)
            {
                redMax = screen[i][j].red;
            }
            if(screen[i][j].green > greenMax)
            {
                greenMax = screen[i][j].green;
            }
            if(screen[i][j].blue > blueMax)
            {
                blueMax = screen[i][j].blue;
            }

            if(screen[i][j].red < redMin)
            {
                redMin = screen[i][j].red;
            }
            if(screen[i][j].green < greenMin)
            {
                greenMin = screen[i][j].green;
            }
            if(screen[i][j].blue < blueMin)
            {
                blueMin = screen[i][j].blue;
            }


        }
    }

    // Apply gamma correction before outputting to file
    double redRange = redMax - redMin;
    double greenRange = greenMax - greenMin;
    double blueRange = blueMax - blueMin;

    double total = (redRange + greenRange + blueRange)/3;

    //std::cout << "red max: " << redMax << " green max: " << greenMax << " blue max: " << blueMax << '\n';
    //std::cout << "red min: " << redMin << " green min: " << greenMin << " blue min: " << blueMin << '\n';

    std::ofstream file;
    file.open(outname.c_str());

    // Writes PPM header
    file << "P3\n128\n128\n255\n";

    // For each column
    for(int j = 0; j<screen_width; j++)
    {
        // Print each value on the row
        for(int i = 0; i<screen_height; i++)
        {
            // Loops through printing RGB values to file
            file << (int) ((screen[i][j].red/total)*255) << " " << (int) ((screen[i][j].green/total)*255) << " " << (int) ((screen[i][j].blue/total)*255) <<" ";
        }

    file << "\n";
    }

    file.close();

}
