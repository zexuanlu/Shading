/*
CSCI 480
Assignment 3 Raytracer

Name: <Your name here>
*/

#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
#include <math.h>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

double pi = 3.14159265;

const int sampling = 2;

double wmin = -(WIDTH * 1.0 / HEIGHT * tan(fov/2 * (pi / 180)));
double hmin = -(tan(fov/2 * (pi/180)));

unsigned char buffer[HEIGHT][WIDTH][3];
unsigned char bufferImage[WIDTH*sampling][HEIGHT*sampling][3];

struct Point {
  double x;
  double y;
  double z;

};

struct Ray {
  struct Point p;
  struct Point d;
};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
bool calRays(Ray &ray, Point &color);
void normalize(Point &unit);

//MODIFY THIS FUNCTION
void draw_scene()
{

  for(int w = 0; w< WIDTH * sampling; w ++){
    for(int h = 0; h < HEIGHT * sampling; h ++){
      Ray ray;
      memset(&ray, 0, sizeof( struct Ray));
      ray.d.x = wmin + 2.0 * w / (WIDTH * sampling * 1.0) * (-wmin);
      ray.d.y = hmin + 2.0 * h / (HEIGHT * sampling * 1.0) * (-hmin);
      ray.d.z = -1;
      normalize(ray.d);
      Point color;
      bool intersection = false;
      intersection = calRays(ray, color);
      if(intersection){ // Put into buffer
        if(color.x > 1)
          bufferImage[w][h][0] = 255;
        else
          bufferImage[w][h][0] = color.x * 255;

        if(color.y > 1)
          bufferImage[w][h][1] = 255;
        else
          bufferImage[w][h][1] = color.y * 255;

        if(color.z > 1)
          bufferImage[w][h][2] = 255;
        else
          bufferImage[w][h][2] = color.z * 255;

      }
      else{
        bufferImage[w][h][0] = 255;
        bufferImage[w][h][1] = 255;
        bufferImage[w][h][2] = 255;
      }
    }
  }

  unsigned int x,y;
  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      int r = 0,g = 0 ,b = 0;
      for(int i = 0; i < sampling; i++){             //supper sampling
        for(int j = 0; j < sampling; j++){
      r += bufferImage[x * sampling + i][y * sampling + j][0];
      g += bufferImage[x * sampling + i][y * sampling + j][1];
      b += bufferImage[x * sampling + i][y * sampling + j][2];
    }
    }
    r = r /pow(sampling, 2);
    g = g /pow(sampling, 2);
    b = b /pow(sampling, 2);

      plot_pixel(x,y,r,g,b);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
  {
    char error[100];
    printf("Expected '%s ' found '%s '\n",expected,found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {

     printf("found triangle\n");
     int j;

     for(j=0;j < 3;j++)
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

void normalize(Point &unit){
  double distance = sqrt(unit.x * unit.x + unit.y * unit.y + unit.z * unit.z);
  unit.x = unit.x / distance;
  unit.y = unit.y / distance;
  unit.z = unit.z / distance;
}

void PointMultiply(Point &a, Point &b, double &result){
  result = a.x * b.x + a.y * b.y + a.z * b.z;

}

void matrixMultiply(double *a, double *b, double *result, int i, int j, int k){
  double s = 0;
  for(int c = 0; c < i; c++){
    for(int d = 0; d< k; d ++){
      for(int e = 0; e < j; e++){
        s = s + a[c * j + e] * b[e * k + d];
      }
      result[c * k + d] = s;
      s = 0;
    }
  }

}

void crossProduct(struct Point * a, struct Point * b, struct Point * result){
  result->x = a->y * b->z - a->z * b->y;
  result->y = a->z * b->x - a->x * b->z;
  result->z = a->x * b->y - a->y * b->x;

}

//Check if the ray intersects with a sphere.
bool intersectCircle(Ray &ray, Sphere &sp, double &t, Point &normal){

  double b = 2*(ray.d.x * (ray.p.x - sp.position[0]) + ray.d.y * (ray.p.y - sp.position[1]) + ray.d.z * (ray.p.z - sp.position[2]));
  double c = (ray.p.x - sp.position[0]) * (ray.p.x - sp.position[0]) + (ray.p.y - sp.position[1]) * (ray.p.y - sp.position[1]) + (ray.p.z - sp.position[2]) * (ray.p.z - sp.position[2]) - sp.radius * sp.radius;
  if((b*b - 4*c) < 1E-8)
    return false; // does not intersect
  double t0 = (-b + sqrt(b*b - 4*c))/2;
  double t1 = (-b - sqrt(b*b - 4*c))/2;
  if(t1 > 1E-8){ //t1 is always smaller, but could be negative
    normal.x = ray.p.x + t1 * ray.d.x - sp.position[0];
    normal.y = ray.p.y + t1 * ray.d.y - sp.position[1];
    normal.z = ray.p.z + t1 * ray.d.z - sp.position[2];
    normalize(normal);
    t = t1;
    return true;
  }
  if(t0 > 1E-8){
    normal.x = ray.p.x + t0 * ray.d.x - sp.position[0];
    normal.y = ray.p.y + t0 * ray.d.y - sp.position[1];
    normal.z = ray.p.z + t0 * ray.d.z - sp.position[2]; // calculate normal
    normalize(normal);
    t = t0;
    return true;

  }

  return false;
}


//Check if the ray intersects with a triangle.
bool intersectTriangle(Ray &ray, Triangle &tr, double &t, double &u, double &v, Point &normal){

  double a[3][3];
  double b[9];
  double c[3];
  double result[3];
  double det;
  a[0][0] = -ray.d.x;
  a[1][0] = -ray.d.y;
  a[2][0] = -ray.d.z;

  a[0][1] = tr.v[1].position[0] - tr.v[0].position[0];
  a[1][1] = tr.v[1].position[1] - tr.v[0].position[1];
  a[2][1] = tr.v[1].position[2] - tr.v[0].position[2];

  a[0][2] = tr.v[2].position[0] - tr.v[0].position[0];
  a[1][2] = tr.v[2].position[1] - tr.v[0].position[1];
  a[2][2] = tr.v[2].position[2] - tr.v[0].position[2];

  det= a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
  b[0]= (a[1][1] * a[2][2] - a[2][1] * a[1][2])/det;
  b[1]= -(a[1][0] * a[2][2] - a[1][2] * a[2][0])/det;
  b[2]= (a[1][0] * a[2][1] - a[2][0] * a[1][1])/det;
  b[3]= -(a[0][1] * a[2][2] - a[0][2] * a[2][1])/det;
  b[4]= (a[0][0] * a[2][2] - a[0][2] * a[2][0])/det;
  b[5]= -(a[0][0] * a[2][1] - a[2][0] * a[0][1])/det;
  b[6]= (a[0][1] * a[1][2] - a[0][2] * a[1][1])/det;
  b[7]= -(a[0][0] * a[1][2] - a[1][0] * a[0][2])/det;
  b[8]= (a[0][0] * a[1][1] - a[1][0] * a[0][1])/det;  // matrix inverse

  c[0] = ray.p.x - tr.v[0].position[0];
  c[1] = ray.p.y - tr.v[0].position[1];
  c[2] = ray.p.z - tr.v[0].position[2];

  matrixMultiply(c, b, result, 1, 3, 3);

  if ( result[0] < 1E-8){
    return false;
  }

  double temp = 1 - result[1] - result[2];
  if(result[1]>= 1E-8 && result[1] <= 1 && result[2] >= 1E-8 && result[2] <= 1 && temp <= 1 && temp >= 1E-8){
    t = result[0];
    u = result[1];
    v = result[2];
    double s = 1 - u - v;
    normal.x = s * tr.v[0].normal[0] + u * tr.v[1].normal[0] + v * tr.v[2].normal[0];
    normal.y = s * tr.v[0].normal[1] + u * tr.v[1].normal[1] + v * tr.v[2].normal[1];
    normal.z = s * tr.v[0].normal[2] + u * tr.v[1].normal[2] + v * tr.v[2].normal[2]; // calculate normal
    return true;
  }
  return false;
}

// Calculate the color.
Point calShading(Ray &ray, Point &normal, double *ks, double *kd, double sh, Point &intersection){
  Point color = {0, 0, 0};
  for( int light_num = 0; light_num < num_lights; light_num++){
    bool inShadow = false;
    Ray shadow;
    shadow.p.x = intersection.x;
    shadow.p.y = intersection.y;
    shadow.p.z = intersection.z;
    shadow.d.x = lights[light_num].position[0] - intersection.x;
    shadow.d.y = lights[light_num].position[1] - intersection.y;
    shadow.d.z = lights[light_num].position[2] - intersection.z;
    normalize(shadow.d);
    Point temp;
    temp.x = lights[light_num].position[0] - intersection.x;
    temp.y = lights[light_num].position[1] - intersection.y;
    temp.z = lights[light_num].position[2] - intersection.z;
    double tempT;
    PointMultiply(temp, temp, tempT);

    for(int i = 0; i < num_spheres; i++){
      double t;
      Point n;
      if(intersectCircle(shadow, spheres[i], t, n)){
        if(t  <= tempT){
          inShadow = true;
        }
      }
    }
                                                                            //shadow rays
    for(int i = 0; i < num_triangles; i++){
      double t, u, v;
      Point n;
      if(intersectTriangle(shadow, triangles[i], t, u, v, n)){
        if(t <= tempT){
          inShadow = true;
        }
      }
    }

    if(!inShadow){
      double dotN;
      PointMultiply(shadow.d, normal, dotN);
      Point r;
      r.x = 2 * dotN * normal.x - shadow.d.x;
      r.y = 2 * dotN * normal.y - shadow.d.y;
      r.z = 2 * dotN * normal.z - shadow.d.z;
      normalize(r);
      if(dotN < 0)
        dotN = 0.0;
      double dotR;
      Point v;
      v.x = - ray.d.x;
      v.y = - ray.d.y;
      v.z = - ray.d.z;
      PointMultiply(r, v, dotR);
      if(dotR < 0)
        dotR = 0.0;
      color.x += lights[light_num].color[0] * (kd[0] * dotN + ks[0] * pow(dotR, sh));
      color.y += lights[light_num].color[1] * (kd[1] * dotN + ks[1] * pow(dotR, sh));
      color.z += lights[light_num].color[2] * (kd[2] * dotN + ks[2] * pow(dotR, sh));
    }
  }
  return color;
}

//Loop through spheres and triangles to calculate the color and shading.
bool calRays(Ray &ray, Point &color){

  bool f = false;
  double distance = 1000000;
  // loops spheres
  for(int i = 0; i < num_spheres; i++){
    double t;
    Point normal;
    if(intersectCircle(ray, spheres[i], t, normal)){
      if(t < distance){
        f = true;
        distance = t;
        Point intersection;
        intersection.x = ray.p.x + t * ray.d.x;
        intersection.y = ray.p.y + t * ray.d.y;
        intersection.z = ray.p.z + t * ray.d.z;
        Point tempColor = calShading(ray, normal, spheres[i].color_specular, spheres[i].color_diffuse, spheres[i].shininess, intersection);

        color.x = tempColor.x;
        color.y = tempColor.y;
        color.z = tempColor.z;

      }

    }

  }
  //loop triangles
  for(int m = 0; m < num_triangles; m++){
    double u, v, t;
    Point normal;

    if(intersectTriangle(ray, triangles[m], t, u, v, normal)){
      if(t < distance){
        distance = t;
        f = true;
        Point intersection;
        intersection.x = ray.p.x + t * ray.d.x;
        intersection.y = ray.p.y + t * ray.d.y;
        intersection.z = ray.p.z + t * ray.d.z;
        double kd[3];
        double ks[3];
        double s = 1 - u - v;
        kd[0] = s * triangles[m].v[0].color_diffuse[0] + u * triangles[m].v[1].color_diffuse[0] + v * triangles[m].v[2].color_diffuse[0];
        kd[1] = s * triangles[m].v[0].color_diffuse[1] + u * triangles[m].v[1].color_diffuse[1] + v * triangles[m].v[2].color_diffuse[1];
        kd[2] = s * triangles[m].v[0].color_diffuse[2] + u * triangles[m].v[1].color_diffuse[2] + v * triangles[m].v[2].color_diffuse[2];
        ks[0] = s * triangles[m].v[0].color_specular[0] + u * triangles[m].v[1].color_specular[0] + v * triangles[m].v[2].color_specular[0];
        ks[1] = s * triangles[m].v[0].color_specular[1] + u * triangles[m].v[1].color_specular[1] + v * triangles[m].v[2].color_specular[1];
        ks[2] = s * triangles[m].v[0].color_specular[2] + u * triangles[m].v[1].color_specular[2] + v * triangles[m].v[2].color_specular[2];
        double sh;
        sh = s * triangles[m].v[0].shininess + u * triangles[m].v[1].shininess + v * triangles[m].v[2].shininess;
        Point temp = calShading(ray, normal, ks, kd, sh, intersection);
        color.x = temp.x;
        color.y = temp.y;
        color.z = temp.z;
      }
    }
  }

  if(f){ // add ambient light at the end
    color.x += ambient_light[0];
    color.y += ambient_light[1];
    color.z += ambient_light[2];
    return true;
  }

  color.x = 0.0;
  color.y = 0.0;
  color.z = 0.0;
  return false;
}




void init()
{
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

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
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
