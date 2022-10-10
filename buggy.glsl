#define MIN_SURFACE_DIST 0.001
#define MAX_DIST 200.
#define PI 3.14159265359
#define ATIME (iTime * 0.2)
#define LOOKAT vec3(-2., 1., -8.)
//#define LOOKFROM ((LOOKAT) + 6.*vec3(cos(ATIME * 2.*PI), 1, sin(-ATIME * 2.*PI)))
#define LOOKFROM vec3(-0., 4., -0.5)

struct Sphere
{
    float r;
};

struct Torus
{
    vec3 p;
    vec2 t;
};

struct Box
{
    vec3 dims;
};

struct Plane
{
    vec3 p;
    vec3 n;
};

float box_sdf(vec3 p, Box c)
{
    vec3 q=abs(p) - c.dims;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float sphere_sdf(vec3 x, in Sphere s)
{
    return length(x) - s.r;
}

float sdTorus(vec3 x, Torus tor)
{
    vec3 p=x - tor.p.xyz;
    vec2 q=vec2(length(p.xz) - tor.t.x, p.y);
    return length(q) - tor.t.y;
}

float sdPlane(in vec3 x, in Plane plane)
{
    return dot(x - plane.p, plane.n);
}

vec2 opU(vec2 a, vec2 b)
{
    return a.x < b.x ? a : b;
}

vec2 opSubtraction(vec2 a, vec2 b)
{
    return -a.x > b.x ? vec2(-a.x, a.y) : b;
}

mat4 translate(vec3 t)
{
    mat4 trans=mat4(1.);
    trans[3]=vec4(-t, 1.);

    return trans;
}

mat4 RotX(float theta)
{
    mat4 trans=mat4(1.);
    trans[1][1]=cos(theta);
    trans[1][2]=sin(theta);
    trans[2][1]=-sin(theta);
    trans[2][2]=cos(theta);

    return transpose(trans);
}

mat4 RotY(float theta)
{
    mat4 trans=mat4(1.);
    trans[0][0]=cos(theta);
    trans[0][2]=-sin(theta);
    trans[2][0]=sin(theta);
    trans[2][2]=cos(theta);

    return transpose(trans);
}

mat4 RotZ(float theta)
{
    mat4 trans=mat4(1.);
    trans[0][0]=cos(theta);
    trans[0][1]=sin(theta);
    trans[1][0]=-sin(theta);
    trans[1][1]=cos(theta);

    return transpose(trans);
}

vec2 get_min_dist(in vec3 p)
{
    float rot_theta=ATIME;
    mat4 rot_x_trans=RotX(rot_theta);
    mat4 rot_y_trans=RotY(rot_theta);
    mat4 rot_z_trans=RotZ(rot_theta * 8.);
    float scale_factor=2.;
    float scale_factor_t=1. / scale_factor;

    Sphere s1=Sphere(.8);
    Sphere s2=Sphere(.8);
    Sphere s3=Sphere(.4);

    Torus tor1=Torus(vec3(0, 0, 0), vec2(1.85, 0.1));
    Box box=Box(vec3(0.2, 0.4, 0.2));
    Plane plane=Plane(vec3(0, -2, 0), normalize(vec3(0, 1, 0.2)));

    //vec2 tmpP = p;
  //  p = 
    mat4 sphere_trans=RotZ(rot_theta * 8.) * translate(vec3(2.5, 2., 0.));

    vec2 tmp=vec2(sphere_sdf((translate(LOOKAT + vec3(0.)) * vec4(p, 1)).xyz, s1), 10.);
    vec2 sdf=opSubtraction(vec2(sphere_sdf((translate(LOOKAT + vec3(0, 0.8, 0.5)) * vec4(p, 1)).xyz, s3), 32.), tmp);

    tmp=vec2(sphere_sdf(p, s2), 51.8);
    sdf=opU(tmp, sdf);
    tmp=vec2(sdPlane(p, plane), 30.);
    sdf=opU(tmp, sdf);
    tmp=vec2(sdTorus((rot_z_trans * translate(vec3(1., 2., -2.)) * vec4(p, 1)).xyz, tor1), 20.);
    sdf=opU(tmp, sdf);
    tmp=vec2(scale_factor * box_sdf(scale_factor_t * (translate(vec3(1., 2., -2)) * vec4(p, 1.)).xyz, box), 40.);
    sdf=opU(tmp, sdf);
    //vec2 sdf= opU(sdf6, opU(sdf5, opU (sdf4, opU (sdf3,  opSubtraction(sdf2, sdf1)))));

    return sdf;
}

bool march(in vec3 ray_origin, in vec3 ray_dir, out vec2 t)
{
    int steps=200;
    float curr=0.;
    for(int i=0;i < steps && curr < MAX_DIST;++i)
    {
        vec3 p=ray_origin + curr * ray_dir;
        vec2 sdf=get_min_dist(p);
        if(abs(sdf.x) < MIN_SURFACE_DIST)
        {
            t=vec2(curr + sdf.x, sdf.y);
            return true;
        }
        curr+=sdf.x;
    }

    return false;
}

vec3 GetNormal(vec3 p)
{
    float d=get_min_dist(p).x;
    vec2 e=vec2(.01, 0);

    vec3 n=d - vec3(get_min_dist(p - e.xyy).x, get_min_dist(p - e.yxy).x, get_min_dist(p - e.yyx).x);

    return normalize(n);
}

vec3 calcNormal(in vec3 p)
{
    const float eps=0.001;
    const vec2 h=vec2(eps, 0);
    return normalize(vec3(get_min_dist(p + h.xyy).x - get_min_dist(p - h.xyy).x, get_min_dist(p + h.yxy).x - get_min_dist(p - h.yxy).x, get_min_dist(p + h.yyx).x - get_min_dist(p - h.yyx).x));
}

vec3 GetLight(in vec3 p, in vec3 n, in float light_comp)
{
    vec3 lookfrom=LOOKFROM;
    vec3 objectColor=0.2 + 0.2 * sin(light_comp * 2. + vec3(0, 1, 2));
    vec3 lightColor=vec3(1.);
    vec3 light_source=vec3(0, 9, 0);

    // ambient
    float ambientStrength=0.1;
    vec3 ambient=ambientStrength * lightColor;

    // diffuse
    vec3 to_light=normalize(light_source - p);
    float diffuse_strength=max(dot(n, to_light), 0.);
    vec3 diffuse=diffuse_strength * lightColor;
    vec2 t;
    if(march(p + n * MIN_SURFACE_DIST * 2., to_light, t))
    {
        float light_distance=length(light_source - p);

        // object is blocking the light source
        if(light_distance > t.x)
        {
            diffuse*=.1;
        }
    }

    // specular
    float specular_brightness=0.5;
    vec3 reflection=normalize(reflect(p - lookfrom, n));
    float specular_value=pow(max(dot(reflection, to_light), 0.), 32.);
    vec3 specular=specular_brightness * specular_value * lightColor;

    // result
    vec3 light=(ambient + diffuse + specular) * objectColor;

    return light;
}

float deg_to_rad(float d)
{
    return (d * PI) / 180.;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    // camera
    float vfov=60.;
    float theta=deg_to_rad(vfov);
    float z_ratio=tan(theta / 2.);
    float ar=iResolution.x / iResolution.y; // ==16:9
    vec3 vup=vec3(0, 1, 0);
    vec3 lookfrom=LOOKFROM;
    vec3 lookat=LOOKAT;
    float th=2. * z_ratio;
    float tw=th * ar;

    vec3 backward=normalize(lookfrom - lookat);
    vec3 right=normalize(cross(vup, backward));
    vec3 up=normalize(cross(backward, right)); // normalize probably unneeded

    vec3 ch=up * th;
    vec3 cw=right * tw;
    vec3 ll=lookfrom - 0.5 * cw - 0.5 * ch - backward;

    float u=(fragCoord.x - .5) / (iResolution.x - 1.);
    float v=(fragCoord.y - .5) / (iResolution.y - 1.);
    vec2 uv=(fragCoord - iResolution.xy * 0.5) / iResolution.y;

    vec3 ray_origin=lookfrom;
    //vec3 ray_dir = normalize(vec3(uv, -1.));
    vec3 ray_dir=normalize(ll + u * cw + v * ch - ray_origin);

    // render

    // background
    fragColor=vec4(0.2, 0.6, 0.5, 1);

    // elements
    vec2 t;
    if(march(ray_origin, ray_dir, t))
    {
        vec3 collision_p=ray_origin + t.x * ray_dir;
        vec3 n=calcNormal(collision_p);
        vec3 col=GetLight(collision_p, n, t.y);// * vec3(0.8, 0.3, 0.2);
        col=pow(col, vec3(.4545));	// gamma correction
        fragColor=vec4(col, 1);
    }
}