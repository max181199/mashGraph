#include <iostream>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include "Bitmap.h"



//// TODO Picture Quadro MF

const uint32_t RED   = 0x000000FF;
const uint32_t GREEN = 0x0000FF00;
const uint32_t BLUE  = 0x00FF0000;
const float kInfinity = std::numeric_limits<float>::max();
int ObjectID = 1;

class Vec3f {
public:
    Vec3f() : x(0), y(0), z(0) {}
    Vec3f(float xx) : x(xx), y(xx), z(xx) {}
    Vec3f(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
    Vec3f operator * (const float &r) const { return Vec3f(x * r, y * r, z * r); }
    Vec3f operator * (const Vec3f &v) const { return Vec3f(x * v.x, y * v.y, z * v.z); }
    Vec3f operator - (const Vec3f &v) const { return Vec3f(x - v.x, y - v.y, z - v.z); }
    Vec3f operator + (const Vec3f &v) const { return Vec3f(x + v.x, y + v.y, z + v.z); }
    Vec3f operator - () const { return Vec3f(-x, -y, -z); }
    Vec3f& operator += (const Vec3f &v) { x += v.x, y += v.y, z += v.z; return *this; }
    friend Vec3f operator * (const float &r, const Vec3f &v)
    { return Vec3f(v.x * r, v.y * r, v.z * r); }
    friend std::ostream & operator << (std::ostream &os, const Vec3f &v)
    { return os << v.x << ", " << v.y << ", " << v.z; }
    float x, y, z;
};
// vector

Vec3f normalize(const Vec3f &v)
{
    float mag2 = v.x * v.x + v.y * v.y + v.z * v.z;
    if (mag2 > 0) {
        float invMag = 1 / sqrtf(mag2);
        return Vec3f(v.x * invMag, v.y * invMag, v.z * invMag);
    }

    return v;
}

inline
float dotProduct(const Vec3f &a, const Vec3f &b)
{ return a.x * b.x + a.y * b.y + a.z * b.z; }
//   (a&b)

Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
    return Vec3f(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
    );
}
//  [a&b]

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0) return false;
    else if (discr == 0) x0 = x1 = - 0.5 * b / a;
    else {
        float q = (b > 0) ?
                  -0.5 * (b + sqrt(discr)) :
                  -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }
    if (x0 > x1) std::swap(x0, x1);
    return true;
}

inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

class Light
{
    public:
        Light(const Vec3f &p, const Vec3f &i) : position(p), intensity(i) {}
        Vec3f position;
        Vec3f intensity; // направление
};

class Eye{
public:
    Eye(const Vec3f &p, const Vec3f &d) : position(p), direction(d) {}
    Vec3f position;
    Vec3f direction;
};

struct Options{
    uint32_t width;
    uint32_t height;
    uint32_t maxDepth = 10;
    uint32_t BackgroundColor  = 0x00000000;
    uint32_t StemRemovalCount = 10;
    uint32_t StemRemovalRadius = 10;
    Vec3f phoneLight = Vec3f(0x00);
    float fov;
    float bias = 1;
};

struct Material{
    float glossFactor = 0.2;
    float transparencyFactor = 0;
    float SpecularFactor = 0;
    float refractionFactor = 0.5;
    Vec3f color = Vec3f(1,1,1);
    float FoneLight = 1;
    float DiffuseLight = 1;
    float SpecularLight = 1;

};


class Object
{
public:
    Object(){
        mymaterial.glossFactor = 300;
        mymaterial.transparencyFactor = 0.2;
        id = ObjectID;
        ObjectID++;
    }
    Object(const Material * mm){
        mymaterial.glossFactor = mm->glossFactor;
        mymaterial.transparencyFactor = mm->transparencyFactor;
        mymaterial.color = mm->color;
        mymaterial.SpecularFactor = mm->SpecularFactor;
        mymaterial.refractionFactor = mm -> refractionFactor;
        mymaterial.DiffuseLight = mm->DiffuseLight;
        mymaterial.SpecularLight = mm->SpecularLight;
        mymaterial.FoneLight = mm->FoneLight;
        id = ObjectID;
        ObjectID++;
    }
    virtual ~Object() {}
    virtual bool intersect(const Vec3f &, const Vec3f &, float &)  = 0;
    virtual Vec3f normal( const Vec3f & Point) = 0;
    int id;
    Material mymaterial;
};




class Sphere : public Object{
public:
    Sphere(const Vec3f &c, const float &r,
           Material * mm) :
       Object(mm) ,center(c), radius(r), radius2(r * r) {}
    bool intersect(const Vec3f & orig, const Vec3f & dir, float & t)
    {
        // analytic solution
        Vec3f L = orig - center;
        float a = dotProduct(dir, dir);
        float b = 2 * dotProduct(dir, L);
        float c = dotProduct(L, L) - radius2;
        float t0, t1;
        if (!solveQuadratic(a, b, c, t0, t1)) {  return false;}
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        t = t0;
        //std :: cout << "Sphere t :: " << t << "\n" << std :: endl;
        return true;
    }

    Vec3f normal( const Vec3f & Point){
        return normalize(Point - center);
    }

    Vec3f center;
    float radius, radius2;
};


class Plain : public Object{
public:
    Plain() :  position( Vec3f(256,0,256) ) , normalP( Vec3f(0,1,0) ) {};
    Plain( const Vec3f &p, const Vec3f &n, Material * mm) :
        Object(mm), position( p ) , normalP( n ) {};
    Plain(const Vec3f &p, const Vec3f &n) : position(p), normalP(n) {};
    Vec3f position;
    Vec3f normalP;
    bool intersect(const Vec3f & orig, const Vec3f & dir, float & t){
        if ( !dotProduct(normalP,dir) ){
            return false;
        }
        double a = dotProduct(dir,normalP);
        Vec3f tmp_vec = orig - position;
        double b = -dotProduct(normalP,tmp_vec);
        float tmp = b/a;
        if( (((dir.z * tmp + orig.z)) >= -1)    && ((dir.z * tmp + orig.z) <= 513)
            && ( (dir.x * tmp + orig.x) >= -1 ) && ( (dir.x * tmp + orig.x) <= 513)
            && ( (dir.y * tmp + orig.y) >= -1)  && ( (dir.y * tmp + orig.y) <= 513)
            ){
            t = tmp;
            return true;
        }
        return false;
    }
    Vec3f normal( const Vec3f & Point){
        return normalize(normalP);
    }

};


class Cube : public Object{
public:
    Cube(const Vec3f & c, float x,float y,float z, Material * mm ) :
            Object(mm) , center(c), x_vector(Vec3f(x,0,0)) ,
            y_vector(Vec3f(0,y,0)), z_vector(Vec3f(0,0,z)) {};

    Vec3f center;
    Vec3f x_vector;
    Vec3f y_vector;
    Vec3f z_vector;
    Vec3f normalL;


    bool intersect(const Vec3f & orig, const Vec3f & dir, float & t)
    {
        float t0;
        float t0_length = kInfinity;
        Vec3f normalP;

        // Part One
        normalP = normalize(crossProduct(z_vector,y_vector));
        if ( dotProduct(crossProduct(y_vector,z_vector),dir) ){
            double a = dotProduct(dir,crossProduct(y_vector,z_vector));
            Vec3f tmp_vec = orig - center;
            double b = -dotProduct(crossProduct(y_vector,z_vector),tmp_vec);
            double tmp = b/a;
            if(
                    ((dir.z * tmp + orig.z) >= center.z) && ((dir.z * tmp + orig.z) <= z_vector.z + center.z) &&
                // ((dir.x * tmp + orig.x) >= 0 ) && ((dir.x * tmp + orig.x) <= 512) &&
                 ((dir.y * tmp + orig.y) >= center.y)  && ( (dir.y * tmp + orig.y) <= y_vector.y + center.y)

                ){
                    if ( dotProduct(orig + dir * tmp,orig + dir * tmp) <= t0_length){
                        t0 = tmp; t0_length =dotProduct(orig + dir * tmp,orig + dir * tmp);
                        normalL = normalP;
                     // normalL = Vec3f(-1,0,0);
                    }
            }
        }
        // Part two
          normalP = normalize(crossProduct(y_vector,x_vector));
        if ( dotProduct(normalP,dir) ){
            double a = dotProduct(dir,normalP);
            Vec3f tmp_vec = orig - center;
            double b = -dotProduct(normalP,tmp_vec);
            double tmp = b/a;
            if(
                    ((dir.z * tmp + orig.z) == center.z) &&
                    ((dir.x * tmp + orig.x) > center.x ) && ((dir.x * tmp + orig.x) < x_vector.x + center.x) &&
                    ((dir.y * tmp + orig.y) > center.y)  && ( (dir.y * tmp + orig.y) < y_vector.y + center.y)

                    ){
                if ( dotProduct(orig + dir * tmp,orig + dir * tmp) < t0_length){
                    t0 = tmp;
                  //  std::cout << "HERE1 = " <<  dotProduct(orig + dir * tmp,orig + dir * tmp) <<  " = " << t0_length << std::endl;
                    t0_length =dotProduct(orig + dir * tmp,orig + dir * tmp);
                    normalL = normalP;
                   normalL = Vec3f(0,0,-1);
                }
            }
        }
//   ////     part3
        normalP = normalize(crossProduct(y_vector,z_vector));
        if ( dotProduct(normalP,dir) ){
            double a = dotProduct(dir,normalP);
            Vec3f tmp_vec = orig - (center+x_vector);
            double b = -dotProduct(normalP,tmp_vec);
            double tmp = b/a;
            if(
                    ((dir.z * tmp + orig.z) >= center.z) && ((dir.z * tmp + orig.z) <= z_vector.z + center.z) &&
                    //((dir.x * tmp + orig.x) >= center.x ) && ((dir.x * tmp + orig.x) <= x_vector.x + center.x) &&
                    ((dir.y * tmp + orig.y) >= center.y)  && ( (dir.y * tmp + orig.y) <= y_vector.y + center.y)

                    ){
                if ( dotProduct(orig + dir * tmp,orig + dir * tmp) < t0_length){
                    t0 = tmp; t0_length =dotProduct(orig + dir * tmp,orig + dir * tmp);
                    normalL = normalP;
                 //   normalL = Vec3f(1,0,0);
                }
            }
        }
////part4
        normalP = normalize(crossProduct(x_vector,z_vector));
        if ( dotProduct(normalP,dir) ){
            double a = dotProduct(dir,normalP);
            Vec3f tmp_vec = orig - center;
            double b = -dotProduct(normalP,tmp_vec);
            double tmp = b/a;
            if(
                ((dir.z * tmp + orig.z) > center.z ) && ((dir.z * tmp + orig.z) < z_vector.z + center.z) &&
                    ((dir.x * tmp + orig.x) > center.x ) && ((dir.x * tmp + orig.x) < x_vector.x + center.x) &&
                        ((dir.y * tmp + orig.y) >= center.y)

                    ){
                if ( dotProduct(orig + dir * tmp,orig + dir * tmp) < t0_length){
                    t0 = tmp; t0_length =dotProduct(orig + dir * tmp,orig + dir * tmp);
                    normalL = normalP;
                    //std::cout << "HERE2 = " <<  dotProduct(orig + dir * tmp,orig + dir * tmp) <<  " = " << t0_length << '\n'<< std::endl;
                    //normalL = Vec3f(0,-1,0);
                }
            }
        }
////        //part5
        normalP = normalize(crossProduct(z_vector,x_vector));
        if ( dotProduct(normalP,dir) ){
            double a = dotProduct(dir,normalP);
            Vec3f tmp_vec = orig - (center + y_vector);
            double b = -dotProduct(normalP,tmp_vec);
            double tmp = b/a;
            if(
                    ((dir.z * tmp + orig.z) >= center.z) && ((dir.z * tmp + orig.z) <= z_vector.z + center.z) &&
                    ((dir.x * tmp + orig.x) >= center.x ) && ((dir.x * tmp + orig.x) <= x_vector.x + center.x) //&&
                // ((dir.y * tmp + orig.y) >= center.y)  && ( (dir.y * tmp + orig.y) <= y_vector.y + center.y)

                    ){
                if ( dotProduct(orig + dir * tmp,orig + dir * tmp) < t0_length){
                    t0 = tmp; t0_length =dotProduct(orig + dir * tmp,orig + dir * tmp);
                    normalL = normalP;
                 //   normalL = Vec3f(0,1,0);
                }
            }
        }
//        //part6
        normalP = normalize(crossProduct(x_vector,y_vector));
        if ( dotProduct(normalP,dir) ){
            double a = dotProduct(dir,normalP);
            Vec3f tmp_vec = orig - (center + z_vector);
            double b = -dotProduct(normalP,tmp_vec);
            double tmp = b/a;
            if(
                    //((dir.z * tmp + orig.z) >= center.z) && ((dir.z * tmp + orig.z) <= z_vector.z + center.z) &&
                    ((dir.x * tmp + orig.x) >= center.x ) && ((dir.x * tmp + orig.x) <= x_vector.x + center.x) &&
                 ((dir.y * tmp + orig.y) >= center.y)  && ( (dir.y * tmp + orig.y) <= y_vector.y + center.y)

                    ){
                if ( (dotProduct(orig + dir * tmp,orig + dir * tmp)) < (t0_length)){
                    t0 = tmp;
                    //std::cout << "HERE2 = " <<  dotProduct(orig + dir * tmp,orig + dir * tmp) << " = " << t0_length << "\n" << std::endl;
                    t0_length =dotProduct(orig + dir * tmp,orig + dir * tmp);
                    normalL=normalP;
                 //   normalL = Vec3f(0,0,1);
                }
            }
        }

        if( t0_length == kInfinity) { return false;}
        else{
            t = t0;
            return true;
        }
    }

    Vec3f normal( const Vec3f & Point){
        return normalL;
    }

};

bool trace(
        const Vec3f &orig, const Vec3f &dir, float & t, uint32_t  index,
        const std::vector<std::unique_ptr<Object>> &objects,Object **hitObject)
{
    *hitObject = nullptr;
    float tmp = kInfinity;
    float t_tmp = t;
    //t = kInfinity;
    //std::cout << t << std::endl;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->intersect(orig, dir ,t_tmp)) {
            Vec3f DeathPoint = orig + dir * t_tmp;
            float Length = dotProduct(DeathPoint,DeathPoint);
            if (
                    (Length < tmp) && (t_tmp >= 0) &&
                            (objects[k]->id != index)
            ){
                //std :: cout << "POINT Z :: " << DeathPoint.z << std :: endl;
                *hitObject = objects[k].get();
                tmp = Length;
                t = t_tmp;
            }
        }
    }
    return (*hitObject != nullptr);
}



Vec3f castRay(
        const Vec3f &orig, const Vec3f &dir,
        const std::vector<std::unique_ptr<Object>> &objects,
        const std::vector<std::unique_ptr<Light>> &lights,
        const Options &options,
        uint32_t index,
        uint32_t depth
        )
{
    if (depth > options.maxDepth) {
        return options.BackgroundColor;
    }

    Vec3f hitColor = Vec3f(options.BackgroundColor);
    Object *hitObject = nullptr;
    float t;
    if (trace(orig, dir,t,index,objects,&hitObject)){
        Vec3f hitPoint = orig + dir * t;
        for (uint32_t i = 0; i < lights.size(); ++i) {
            Vec3f lightDir = lights[i]->position - hitPoint;
            Vec3f lightIntencity = Vec3f( lights[i]->intensity.x,lights[i]->intensity.y,lights[i]->intensity.z);
            Vec3f NN = hitObject->normal(hitPoint);
            Vec3f pix = Vec3f(0,0,0);
            Vec3f refract = Vec3f(0,0,0);

            //REFRACTION BLOCK
            if (hitObject->mymaterial.transparencyFactor){
                float tmp = 1 ;
                if (hitObject->id == (-index)){
                    tmp = 1/hitObject->mymaterial.refractionFactor;
                } else{
                    tmp = hitObject->mymaterial.refractionFactor;
                }
                float L = float(hitObject->mymaterial.refractionFactor*tmp)
                        - float(1);
                L = L / ( dotProduct(dir,NN) * dotProduct(dir,NN) );
                L = L+1;
                L = sqrt(L) - 1;
                L = L * dotProduct(dir,NN);
                Vec3f RefractionDirection = normalize( dir + NN * L);
                Vec3f RefractionOrigin    = hitPoint + RefractionDirection*0.01;
                if ( hitObject->id == (-index) ){
                    refract = castRay(RefractionOrigin,RefractionDirection,objects,lights,options,0,depth+1);
                } else {
                    refract = castRay(RefractionOrigin,RefractionDirection,objects,lights,options,-hitObject->id,depth+1);
                }
            }
            //REFRACTION BLOCK
            //SPECULAR BLOCK
            if (hitObject->mymaterial.SpecularFactor){
                Vec3f SpecularDirection = normalize(2 * NN + dir);
                Vec3f SpecularOrigin    = hitPoint ;
                pix = castRay(SpecularOrigin,SpecularDirection,objects,lights,options,hitObject->id,depth + 1);
            }
            //SPECULAR BLOCK
            // SHADOW BLOCK
                Object *shadowHitObject = nullptr;
                for (uint32_t k = 0; k < objects.size(); ++k) {
                    if (objects[k]->intersect(hitPoint, lightDir ,t)) {
                        shadowHitObject = objects[k].get();
                        Vec3f shadowHitPoint  = hitPoint + lightDir * t;
                        if ( shadowHitObject == hitObject ) { continue;}
                        if ( dotProduct(lightDir,lightDir) <= dotProduct(lights[i]->position - shadowHitPoint
                                ,lights[i]->position - shadowHitPoint)) { continue;}
                        if ( dotProduct(lightDir, normalize(shadowHitPoint - lights[i]->position)) > 0 ){ continue;}
                        lightIntencity.x = lightIntencity.x * shadowHitObject->mymaterial.transparencyFactor;
                        lightIntencity.y = lightIntencity.y * shadowHitObject->mymaterial.transparencyFactor;
                        lightIntencity.z = lightIntencity.z * shadowHitObject->mymaterial.transparencyFactor;
                        pix.x = pix.x * shadowHitObject->mymaterial.transparencyFactor;
                        pix.y = pix.y * shadowHitObject->mymaterial.transparencyFactor;
                        pix.z = pix.z * shadowHitObject->mymaterial.transparencyFactor;
                        if(shadowHitObject->mymaterial.transparencyFactor == 0) break;
                   }
                }
            // SHADOW BLOCK
            // TRANSPARENCY BLOCK
//            Object *TranHitObject = nullptr;
//                for (uint32_t k = 0; k < objects.size(); ++k) {
//                    if (objects[k]->intersect(hitPoint, lightDir ,t)) {
//                        shadowHitObject = objects[k].get();
//                        Vec3f shadowHitPoint  = hitPoint + lightDir * t;
//                        if ( shadowHitObject == hitObject ) { continue;}
//                        if ( dotProduct(lightDir,lightDir) <= dotProduct(lights[i]->position - shadowHitPoint
//                                ,lights[i]->position - shadowHitPoint)) { continue;}
//                        if ( dotProduct(lightDir, normalize(shadowHitPoint - lights[i]->position)) > 0 ){ continue;}
//                        lightIntencity.x = lightIntencity.x * shadowHitObject->mymaterial.transparencyFactor;
//                        lightIntencity.y = lightIntencity.y * shadowHitObject->mymaterial.transparencyFactor;
//                        lightIntencity.z = lightIntencity.z * shadowHitObject->mymaterial.transparencyFactor;
//                        pix.x = pix.x * shadowHitObject->mymaterial.transparencyFactor;
//                        pix.y = pix.y * shadowHitObject->mymaterial.transparencyFactor;
//                        pix.z = pix.z * shadowHitObject->mymaterial.transparencyFactor;
//                        if(shadowHitObject->mymaterial.transparencyFactor == 0) break;
//                   }
//                }

            // TRANSPARENCY BLOCK
            float lightDistance2 = sqrt(dotProduct(lightDir, lightDir)) ;
            lightDistance2 = (sqrt(lightDistance2));
            lightDir = normalize(lightDir);

            // Фоновое освещение
            hitColor += options.phoneLight * (hitObject->mymaterial.FoneLight);
            // Фоновое освещение

            // Рассеянный свет
            float Lambert = dotProduct(lightDir,NN);
            if (Lambert < 0) {
                Lambert = 0;
            }
            hitColor += Vec3f(  (Lambert/lightDistance2*lightIntencity.x * hitObject->mymaterial.color.x * hitObject->mymaterial.DiffuseLight),
                    (Lambert/lightDistance2*lightIntencity.y* hitObject->mymaterial.color.y * hitObject->mymaterial.DiffuseLight),
                    (Lambert/lightDistance2*lightIntencity.z) * hitObject->mymaterial.color.z * hitObject->mymaterial.DiffuseLight);
            // Рассеянный свет

            // Зеркальный свет
            float mirrorLight =pow( dotProduct(NN,normalize(  lightDir+ dir * -1)),hitObject->mymaterial.glossFactor) * hitObject->mymaterial.SpecularLight;
            if (mirrorLight < 0) { mirrorLight = 0;}
            hitColor += Vec3f(  (mirrorLight/lightDistance2*lightIntencity.x* hitObject->mymaterial.color.x),
                               (mirrorLight/lightDistance2*lightIntencity.y* hitObject->mymaterial.color.y),
                               (mirrorLight/lightDistance2*lightIntencity.z)* hitObject->mymaterial.color.z );
            // Зеркальный свет


            if ( (hitObject->mymaterial.SpecularFactor != 0) || (hitObject->mymaterial.transparencyFactor !=0) ){
               hitColor =hitColor * ( 1 - hitObject->mymaterial.SpecularFactor - hitObject->mymaterial.transparencyFactor) +
                       pix * hitObject->mymaterial.SpecularFactor +
                       refract * hitObject->mymaterial.transparencyFactor ;
            }

            if (hitColor.x > 0xFF) { hitColor.x = 0xFF;}
            if (hitColor.y > 0xFF) { hitColor.y = 0xFF;}
            if (hitColor.z > 0xFF) { hitColor.z = 0xFF;}
            if (hitColor.x < 0x00) { hitColor.x = 0x00;}
            if (hitColor.y < 0x00) { hitColor.y = 0x00;}
            if (hitColor.z < 0x00) { hitColor.z = 0x00;}
        }
    }
    return hitColor;
}

//std::cout<< " hitColor(x: " << hitColor.x << ", y: "<< hitColor.y << ", z: " << hitColor.z << ")\n" << std::endl; // TEST

void render(
        std::vector<uint32_t> &image,
        const Options &options,
        const std::vector<std::unique_ptr<Object>> &objects,
        const std::vector<std::unique_ptr<Light>> &lights)
{
     // Our Picture
    for(auto& pixel : image)
        pixel = 0x00000000; // Create Black Phone
    Vec3f pix(0xFF);
    Vec3f orig(256,256,-100);
    for (uint32_t j = 0; j < options.height; ++j) {
        for (uint32_t i = 0; i < options.width; ++i) {
            // generate primary ray direction
            float x = ( float(i) - options.width/2)   ;
            float y = ( float(j) - options.height/2)   ;
            Vec3f dir = normalize(Vec3f(x, y, 255));
            uint32_t index = -100;
            pix = castRay(orig, dir, objects, lights, options,index, 1);

            // Устранение ступенчатости
            Vec3f Step = pix;
                for ( int m = 0; m < 2 * options.StemRemovalCount; m++ ){
                    for ( int n = 0; n < 2 * options.StemRemovalCount; n++ ){
                        int Xx = m - options.StemRemovalCount;
                        int Yy = n - options.StemRemovalCount ;
                        float XXx = Xx;
                        float YYy = Yy;
                        XXx=XXx/options.StemRemovalRadius;
                        YYy=YYy/options.StemRemovalRadius;
                        dir = normalize(Vec3f(x+ XXx,y+YYy,255));
                        Vec3f one = Vec3f(255);
                        one = castRay(orig, dir, objects, lights, options,index, 1);
                        Step += one;
                    }
                }
                float tmp = 1;
                if ( options.StemRemovalCount){
                    tmp = tmp / (options.StemRemovalCount * options.StemRemovalCount * 4 + 1);
                }
                Step = Step * tmp;
                pix = Step;
            // Устранение ступенчатости

            image[ j*options.height + i  ] = (int(pix.y) << 8) + (int(pix.x) << 16) + int(pix.z) ;
        }
    }
}




int main(int argc, const char** argv)
{
  std::unordered_map<std::string, std::string> cmdLineParams;

  for(int i=0; i<argc; i++)
  {
    std::string key(argv[i]);

    if(key.size() > 0 && key[0]=='-')
    {
      if(i != argc-1) // not last argument
      {
        cmdLineParams[key] = argv[i+1];
        i++;
      }
      else
        cmdLineParams[key] = "";
    }
  }

  std::string outFilePath = "zout.bmp";
  if(cmdLineParams.find("-out") != cmdLineParams.end())
    outFilePath = cmdLineParams["-out"];

  int sceneId = 0;
  if(cmdLineParams.find("-scene") != cmdLineParams.end())
    sceneId = atoi(cmdLineParams["-scene"].c_str());

  Material Basic;
  Basic.glossFactor = 10;
  Basic.color = Vec3f(1,1,1);
  Basic.transparencyFactor = 0;
  Basic.refractionFactor = 0;
  Basic.SpecularFactor = 0;
  Basic.DiffuseLight = 1;
  Basic.FoneLight = 1;
  Basic.SpecularLight = 0.1;

  Material Mirror;
  Mirror.glossFactor = 10;
  Mirror.color = Vec3f(1,1,1);
  Mirror.transparencyFactor = 0;
  Mirror.refractionFactor = 0;
  Mirror.SpecularFactor = 1;
  Mirror.DiffuseLight = 0;
  Mirror.FoneLight = 0;
  Mirror.SpecularLight = 0;

  Material UnrealGlass;
  UnrealGlass.glossFactor=10;
  UnrealGlass.color = Vec3f (1,1,1);
  UnrealGlass.transparencyFactor = 0.9;
  UnrealGlass.refractionFactor = 1.5;
  UnrealGlass.SpecularFactor = 0;
  UnrealGlass.DiffuseLight = 1;
  UnrealGlass.FoneLight = 1;
  UnrealGlass.SpecularLight = 0.1;

    std::vector<std::unique_ptr<Object>> objects;
    std::vector<std::unique_ptr<Light>> lights;

    Plain * plain = new  Plain( Vec3f(256,0,256),Vec3f(0,1,0),&Basic);
    Basic.color = Vec3f(1,0.7,0.4);
    Plain * plaint_test = new Plain( Vec3f(256,256,510),Vec3f(0,0,-1),&Basic);
    Plain * plaint_test2 = new Plain ( Vec3f(0,256,256), Vec3f(1,0,0),&Mirror);
    Basic.color = Vec3f(0.7,1,0.4);
    Plain * plaint_test3 = new Plain ( Vec3f(512,256,256), Vec3f(-1,0,0),&Basic);
    Basic.color = Vec3f(1,1,1);
    Plain * plaint_test4 = new Plain ( Vec3f(256,512,256), Vec3f(0,-1,0),&Basic);
    Basic.color = Vec3f ( 1,1,1);
    Cube  * cube_test    = new Cube(Vec3f(150,100,300),70,70,70,&Basic);
    Basic.color = Vec3f ( 1,1,1);
    Sphere * sphere = new Sphere(Vec3f(400,110,400),70,&UnrealGlass);
    Sphere * sphere2 = new Sphere(Vec3f(481,170,400),30,&UnrealGlass);



    objects.push_back(std::unique_ptr<Plain>(plaint_test));
    objects.push_back(std::unique_ptr<Plain>(plaint_test2));
    objects.push_back(std::unique_ptr<Plain>(plaint_test3));
    objects.push_back(std::unique_ptr<Plain>(plaint_test4));
    objects.push_back(std::unique_ptr<Plain>(plain));
    objects.push_back(std::unique_ptr<Cube>(cube_test));
    objects.push_back(std::unique_ptr<Sphere>(sphere));
    objects.push_back(std::unique_ptr<Sphere>(sphere2));

    lights.push_back(std::unique_ptr<Light>(new Light(Vec3f(256, 400, 256),
            Vec3f(3000,3000,3000)                )));
//    lights.push_back(std::unique_ptr<Light>(new Light(Vec3f(256, 400, 100),
//                                                      Vec3f(3000,0,0)                )));
//
//    lights.push_back(std::unique_ptr<Light>(new Light(Vec3f(100, 300, 100),
//                                                      Vec3f(0,3000,0)                )));
//
//    lights.push_back(std::unique_ptr<Light>(new Light(Vec3f(400, 50, 100),
//                                                      Vec3f(0,0,3000)                )));


    Options options;
    options.width = 512;
    options.height = 512;
    options.BackgroundColor = 0x00;
    options.phoneLight = Vec3f(0x12);
    options.maxDepth = 10;
    options.StemRemovalCount = 1;
    options.StemRemovalRadius = 10;

    std::vector<uint32_t> image(options.height * options.width);

    render(image,options, objects, lights);

    SaveBMP(outFilePath.c_str(), image.data(), options.width, options.height);



    std::cout << "end." << std::endl;
  return 0;
}