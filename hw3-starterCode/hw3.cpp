/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Aishwarya Verghese
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <glm/glm.hpp>
#include <imageIO.h>

#include <algorithm>

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

using namespace std;

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY; //default - no image

#define WIDTH 640
#define HEIGHT 480


#define fov 60.0 //the field of view of the camera
const double PI = 3.141592f;
double aspect_ratio = (double)WIDTH/(double)HEIGHT;
double angle = std::tan ((fov / 2.0) * (PI / 180.0));

//Calculating image plane coordinates from feild of view and aspect ratio
double left_x = - aspect_ratio * angle;
double bottom_y = -angle;

double x_one_step = (2 * aspect_ratio * angle) / (double)WIDTH;
double y_one_step = (2 * angle) / (double)HEIGHT;

//Buffer to store image data
unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

glm::dvec3 camera_ray_origin = glm::dvec3(0.0,0.0,0.0); //Assuming camera at origin

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//Get all triangle intersection vertex properties - position, normal, diffused color, specular color and shininess
Vertex getTriangleIntersectionVertex(float alpha,float beta,float gamma,Triangle triangle)
{

  Vertex intersection;

  intersection.position[0] = (alpha * triangle.v[0].position[0]) + (beta * triangle.v[1].position[0]) + (gamma * triangle.v[2].position[0]);
  intersection.position[1] = (alpha * triangle.v[0].position[1]) + (beta * triangle.v[1].position[1]) + (gamma * triangle.v[2].position[1]);
  intersection.position[2] = (alpha * triangle.v[0].position[2]) + (beta * triangle.v[1].position[2]) + (gamma * triangle.v[2].position[2]);

  intersection.normal[0] = (alpha * triangle.v[0].normal[0]) + (beta * triangle.v[1].normal[0]) + (gamma * triangle.v[2].normal[0]);
  intersection.normal[1] = (alpha * triangle.v[0].normal[1]) + (beta * triangle.v[1].normal[1]) + (gamma * triangle.v[2].normal[1]);
  intersection.normal[2] = (alpha * triangle.v[0].normal[2]) + (beta * triangle.v[1].normal[2]) + (gamma * triangle.v[2].normal[2]);
  glm::dvec3 normalized_normal = glm::normalize(glm::dvec3(intersection.normal[0],intersection.normal[1],intersection.normal[2]));
  intersection.normal[0] = normalized_normal.x;
  intersection.normal[1] = normalized_normal.y;
  intersection.normal[2] = normalized_normal.z;

  intersection.color_diffuse[0] = (alpha * triangle.v[0].color_diffuse[0]) + (beta * triangle.v[1].color_diffuse[0]) + (gamma * triangle.v[2].color_diffuse[0]);
  intersection.color_diffuse[1] = (alpha * triangle.v[0].color_diffuse[1]) + (beta * triangle.v[1].color_diffuse[1]) + (gamma * triangle.v[2].color_diffuse[1]);
  intersection.color_diffuse[2] = (alpha * triangle.v[0].color_diffuse[2]) + (beta * triangle.v[1].color_diffuse[2]) + (gamma * triangle.v[2].color_diffuse[2]);

  intersection.color_specular[0] = (alpha * triangle.v[0].color_specular[0]) + (beta * triangle.v[1].color_specular[0]) + (gamma * triangle.v[2].color_specular[0]);
  intersection.color_specular[1] = (alpha * triangle.v[0].color_specular[1]) + (beta * triangle.v[1].color_specular[1]) + (gamma * triangle.v[2].color_specular[1]);
  intersection.color_specular[2] = (alpha * triangle.v[0].color_specular[2]) + (beta * triangle.v[1].color_specular[2]) + (gamma * triangle.v[2].color_specular[2]);

  intersection.shininess = (alpha * triangle.v[0].shininess) + (beta * triangle.v[1].shininess) + (gamma * triangle.v[2].shininess);

  return intersection;
}

//Get all sphere intersection vertex properties - position, normal, diffused color, specular color and shininess
Vertex getSphereIntersectionVertex(Sphere sphere, glm::dvec3 ray_origin, glm::dvec3 ray_direction, double t)
{
  Vertex intersection;

  intersection.color_diffuse[0] = sphere.color_diffuse[0];
  intersection.color_diffuse[1] = sphere.color_diffuse[1];
  intersection.color_diffuse[2] = sphere.color_diffuse[2];
  intersection.color_specular[0] = sphere.color_specular[0];
  intersection.color_specular[1] = sphere.color_specular[1];
  intersection.color_specular[2] = sphere.color_specular[2];
  intersection.shininess = sphere.shininess;

  intersection.position[0] = ray_origin.x + ray_direction.x * t;
  intersection.position[1] = ray_origin.y + ray_direction.y * t;
  intersection.position[2] = ray_origin.z + ray_direction.z * t;

  intersection.normal[0] = (intersection.position[0] - sphere.position[0]) / sphere.radius;
  intersection.normal[1] = (intersection.position[1] - sphere.position[1]) / sphere.radius;
  intersection.normal[2] = (intersection.position[2] - sphere.position[2]) / sphere.radius;

  return intersection;
}

//Get Phong lighting color given one light source and point
glm::dvec3 getPhongLighting(Light light,Vertex point)
{
  //Calculate light and view direction vectors
  glm::dvec3 light_dir = glm::normalize(glm::dvec3((light.position[0] - point.position[0]),(light.position[1] - point.position[1]),(light.position[2] - point.position[2])));
  glm::dvec3 view = glm::normalize(glm::dvec3(0 - point.position[0],0 - point.position[1],0 - point.position[2]));

  //Calculated reflected vector
  glm::dvec3 reflect_dir = -glm::reflect(light_dir, glm::dvec3(point.normal[0],point.normal[1],point.normal[2]));
  reflect_dir = glm::normalize(reflect_dir);

  //Phong Lighting calculations
  double kd = max(glm::dot(light_dir,glm::dvec3(point.normal[0],point.normal[1],point.normal[2])) , 0.0);
  double ks = max(glm::dot(reflect_dir,view) , 0.0);
  double r = light.color[0] * ((point.color_diffuse[0] * kd) +  (point.color_specular[0] * pow(ks,point.shininess)));
  double g = light.color[1] * ((point.color_diffuse[1] * kd) +  (point.color_specular[1] * pow(ks,point.shininess)));
  double b = light.color[2] * ((point.color_diffuse[2] * kd) +  (point.color_specular[2] * pow(ks,point.shininess)));

  return glm::dvec3(r,g,b);
}

//Checks if ray intersects triangle
bool triangleIntersect(glm::dvec3 normal,glm::dvec3 a,glm::dvec3 b,glm::dvec3 c,glm::dvec3 ray_origin, glm::dvec3 ray_direction,float &alpha,float &beta,float &gamma,double &depth)
{
  //Check if ray is parallel to normal
  double parallel = glm::dot(normal,ray_direction);
  if(parallel == 0) //ray parallel to triangle
    return false;


  double d = - glm::dot(normal,a);
  double t = - (glm::dot(normal,ray_origin) + d) / parallel;
  if(t <= 0)  //Doesn't intersect plane
  {
    return false;
  }

  //get intersection point of ray and plane containing triangle
  glm::dvec3 i;
  i.x = ray_origin.x + ray_direction.x * t;
  i.y = ray_origin.y + ray_direction.y * t;
  i.z = ray_origin.z + ray_direction.z * t;
  depth = i.z;

  //check if point in trangle
  float abc_area, abi_area, aic_area, ibc_area;
  if(glm::dot(normal,glm::dvec3(0,0,1.0)) != 0) // if XY plane not perpendicular to triangle
  { //Project to XY plane
    abc_area = 0.5 * (((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y)));
    abi_area = 0.5 * (((b.x - a.x) * (i.y - a.y)) - ((i.x - a.x) * (b.y - a.y)));
    aic_area = 0.5 * (((i.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (i.y - a.y)));
    ibc_area = 0.5 * (((b.x - i.x) * (c.y - i.y)) - ((c.x - i.x) * (b.y - i.y)));
  }
  else if(glm::dot(normal,glm::dvec3(0,1.0,0)) != 0)//XZ plane not perpendicular to triangle
  {
    //Project to XZ plane
    abc_area = 0.5 * (((b.x - a.x) * (c.z - a.z)) - ((c.x - a.x) * (b.z - a.z)));
    abi_area = 0.5 * (((b.x - a.x) * (i.z - a.z)) - ((i.x - a.x) * (b.z - a.z)));
    aic_area = 0.5 * (((i.x - a.x) * (c.z - a.z)) - ((c.x - a.x) * (i.z - a.z)));
    ibc_area = 0.5 * (((b.x - i.x) * (c.z - i.z)) - ((c.x - i.x) * (b.z - i.z)));
  }
  else
  {
    //Project to YZ plane
    abc_area = 0.5 * (((b.y - a.y) * (c.z - a.z)) - ((c.y - a.y) * (b.z - a.z)));
    abi_area = 0.5 * (((b.y - a.y) * (i.z - a.z)) - ((i.y - a.y) * (b.z - a.z)));
    aic_area = 0.5 * (((i.y - a.y) * (c.z - a.z)) - ((c.y - a.y) * (i.z - a.z)));
    ibc_area = 0.5 * (((b.y - i.y) * (c.z - i.z)) - ((c.y - i.y) * (b.z - i.z)));
  }

  alpha = ibc_area / abc_area;
  beta = aic_area / abc_area;
  gamma = abi_area / abc_area;

  //if point lies in triangle
  if( (alpha >= 0.0) && (alpha <= 1.0)  && (beta >= 0.0) && (beta <= 1.0) && (gamma >= 0.0) && (gamma <= 1.0))
    return true;
  else
    return false;
}

//Checks if ray intersects sphere
bool sphereIntersect(Sphere sphere, glm::dvec3 ray_origin, glm::dvec3 ray_direction, double &t)
{
  double b = 2 * (ray_direction.x * (ray_origin.x - sphere.position[0]) + ray_direction.y * (ray_origin.y - sphere.position[1]) + ray_direction.z * (ray_origin.z - sphere.position[2]));
  double c = pow((ray_origin.x - sphere.position[0]),2) + pow((ray_origin.y - sphere.position[1]),2) + pow((ray_origin.z - sphere.position[2]),2) - sphere.radius*sphere.radius;
  double disc = b*b - 4*c;  //calculate discrimnant
  if(disc == 0.0) //intersects 1 point
  {
    t = -b/2;
    return true;
  }
  else if(disc > 0.0) //instersects 2 points
  {
    double t0 = -(b + sqrt(disc)) / 2.0;
    double t1 = -(b - sqrt(disc)) / 2.0;
    if(t0 <= 0.0 || t1 <= 0.0) //intersects backwards
      return false;

    t = min(t0,t1); //get closest of the two
    return true;
  }
  return false;
}

//returns true if point is in shadow with respect to a light
bool inShadow(Light light,Vertex point,int this_triangle_num,int this_sphere_num)
{
    //Calculate shadow ray vector
    glm::dvec3 shadow_ray = glm::normalize(glm::dvec3(light.position[0] - point.position[0], light.position[1] - point.position[1], light.position[2] - point.position[2]));

    //for all triangles
    for(int triangle_num = 0 ; triangle_num < num_triangles; triangle_num++)
    {
      if(this_triangle_num == triangle_num) //exclude object firing the shadow ray
        continue;

      float alpha,beta,gamma;
      double depth;
      //Get all vertices and normal
      glm::dvec3 a = glm::dvec3(triangles[triangle_num].v[0].position[0],triangles[triangle_num].v[0].position[1],triangles[triangle_num].v[0].position[2]);
      glm::dvec3 b =  glm::dvec3(triangles[triangle_num].v[1].position[0],triangles[triangle_num].v[1].position[1],triangles[triangle_num].v[1].position[2]);
      glm::dvec3 c =  glm::dvec3(triangles[triangle_num].v[2].position[0],triangles[triangle_num].v[2].position[1],triangles[triangle_num].v[2].position[2]);
      glm::dvec3 normal = glm::normalize(glm::cross((b-a),(c-a)));

      //Does shadow ray intersect the triangle?
      if(triangleIntersect(normal,a,b,c,glm::dvec3(point.position[0],point.position[1],point.position[2]),shadow_ray,alpha,beta,gamma,depth))
      {
        Vertex intersection;
        intersection.position[0] = (alpha * triangles[triangle_num].v[0].position[0]) + (beta * triangles[triangle_num].v[1].position[0]) + (gamma * triangles[triangle_num].v[2].position[0]);
        intersection.position[1] = (alpha * triangles[triangle_num].v[0].position[1]) + (beta * triangles[triangle_num].v[1].position[1]) + (gamma * triangles[triangle_num].v[2].position[1]);
        intersection.position[2] = (alpha * triangles[triangle_num].v[0].position[2]) + (beta * triangles[triangle_num].v[1].position[2]) + (gamma * triangles[triangle_num].v[2].position[2]);

        //If the distance between intersecting object and source object is shorter than the distance between the light source and the source object then the source object is in shadow
        double light_dist = glm::distance(glm::dvec3(point.position[0],point.position[1],point.position[2]),glm::dvec3(light.position[0],light.position[1],light.position[2]));
        double obj_dist = glm::distance(glm::dvec3(point.position[0],point.position[1],point.position[2]),glm::dvec3(intersection.position[0],intersection.position[1],intersection.position[2]));

        if(light_dist > obj_dist)
          return true;
      }
    }

    //for all the spheres
    for(int sphere_num = 0; sphere_num < num_spheres; sphere_num++)
    {
      if(this_sphere_num == sphere_num) //exclude source sphere
        continue;

      Sphere sphere = spheres[sphere_num];
      double t;
      if(sphereIntersect(sphere, glm::dvec3(point.position[0],point.position[1],point.position[2]), shadow_ray, t))
      {
        double sphere_i[3];
        sphere_i[0] = point.position[0] + shadow_ray.x * t;
        sphere_i[1] = point.position[1] + shadow_ray.y * t;
        sphere_i[2] = point.position[2] + shadow_ray.z * t;

        //If the distance between intersecting object and source object is shorter than the distance between the light source and the source object then the source object is in shadow
        double light_dist = glm::distance(glm::dvec3(point.position[0],point.position[1],point.position[2]),glm::dvec3(light.position[0],light.position[1],light.position[2]));
        double obj_dist = glm::distance(glm::dvec3(point.position[0],point.position[1],point.position[2]),glm::dvec3(sphere_i[0],sphere_i[1],sphere_i[2]));

        if(light_dist > obj_dist)
          return true;
      }
    }
    return false; //No intersecting objects so not in shadow
}

void draw_scene()
{

  //For each row
  double x_pixel = left_x;
  for(unsigned int x = 0 ; x < WIDTH; x++)
  {
    glPointSize(2.0);
    glBegin(GL_POINTS);

    //For each column
    double y_pixel = bottom_y;
    for(unsigned int y = 0 ; y < HEIGHT; y++)
    {
      //initalizations
      glm::dvec3 ray_direction = glm::normalize(glm::dvec3(x_pixel,y_pixel,-1.0));
      glm::dvec3 final_color = glm::dvec3(0.0,0.0,0.0);
      double closest_z = +1.0;
      int closest_triangle = -1;
      int closest_sphere = -1;
      Vertex final_intersection;

      //For each triangle
      for(int triangle_num = 0; triangle_num < num_triangles; triangle_num++)
      {
        float alpha,beta,gamma;
        double depth;
        //Get all vertices and normal
        glm::dvec3 a = glm::dvec3(triangles[triangle_num].v[0].position[0],triangles[triangle_num].v[0].position[1],triangles[triangle_num].v[0].position[2]);
        glm::dvec3 b =  glm::dvec3(triangles[triangle_num].v[1].position[0],triangles[triangle_num].v[1].position[1],triangles[triangle_num].v[1].position[2]);
        glm::dvec3 c =  glm::dvec3(triangles[triangle_num].v[2].position[0],triangles[triangle_num].v[2].position[1],triangles[triangle_num].v[2].position[2]);
        glm::dvec3 normal = glm::normalize(glm::cross((b-a),(c-a)));

         //Does ray intersect the triangle?
        if(!triangleIntersect(normal,a,b,c,camera_ray_origin,ray_direction,alpha,beta,gamma,depth))
          continue;
        else
        {
          //calculate if closest intersection
          if(closest_z == 1.0 || depth > closest_z)
          {
            final_intersection = getTriangleIntersectionVertex(alpha,beta,gamma,triangles[triangle_num]);
            closest_z = depth;  //update closest z value and object number
            closest_triangle = triangle_num;
            closest_sphere = -1;
          }
        }
      }

      //do it for each sphere
      for(int sphere_num = 0; sphere_num < num_spheres; sphere_num++)
      {
        Sphere sphere  = spheres[sphere_num];
        double t;
        if(sphereIntersect(sphere, camera_ray_origin, ray_direction, t))
        {
          double depth = camera_ray_origin.z + ray_direction.z * t;

          if(closest_z == 1.0 || depth > closest_z)
          {
            final_intersection = getSphereIntersectionVertex(sphere, camera_ray_origin, ray_direction, t);
            closest_z = depth;
            closest_triangle = -1;
            closest_sphere = sphere_num;
          }
        }
      }

      if(closest_z != 1.0)  //intersected some object
      {
        //for each light source
        for(int light_num = 0; light_num <= num_lights;light_num++)
        {
          if(!inShadow(lights[light_num],final_intersection,closest_triangle,closest_sphere))
          {
            glm::dvec3 color = getPhongLighting(lights[light_num],final_intersection);
            final_color.r += color.r;
            final_color.b += color.b;
            final_color.g += color.g;
          }
        }
        //add ambient lights & clip
        final_color.r = std::min((final_color.r + ambient_light[0]) , 1.0);
        final_color.g = std::min((final_color.g + ambient_light[1]) , 1.0);
        final_color.b = std::min((final_color.b + ambient_light[2]) , 1.0);
        plot_pixel(x,y,final_color.r * 255.0,final_color.g * 255.0,final_color.b * 255.0);
      }
      else
        plot_pixel(x,y,255.0,255.0,255.0);

      y_pixel = y_pixel + y_one_step;
    }
    glEnd();
    glFlush();
    x_pixel = x_pixel + x_one_step;
  }
  printf("Done!\n"); fflush(stdout);
}

//Plot pixel to display
void plot_pixel_display(int x, int y, unsigned char r,unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

//Plot pixel to image (buffer)
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

//After calculations save image buffer to JPEG file
void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else
    printf("File saved Successfully\n");
}

//Parse values from input file
void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

//Load values from input file
int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  //Initial transformation matrices
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3) //draw to image and display
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2) //draw to display only
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
