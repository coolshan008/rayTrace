#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <time.h>
#define PI 3.141592
#define HEIGHT 300
#define WEIGHT 300
#define NUMBEROFMODELS 3

//int Count = 0;
//int Count2 = 0;
/*
屏幕上下为y（上正），左右为x（右正），往外为z正*/

using namespace std;


//fstream out_file("debug.txt", ios::out);

struct Point{
	GLfloat x, y, z;
	GLfloat normx, normy, normz;
};

struct Triangle{
	glm::vec3 P1, P2, P3;
	glm::vec3 norm;
};

struct Index{
	GLuint xindex, yindex, zindex;
};

struct interPointIndex{
	GLint model, triangle;
	GLfloat distance;
};

bool IntersectTriangle(glm::vec3 &orig, glm::vec3 &dir,
	glm::vec3 &v0, glm::vec3 &v1, glm::vec3 &v2,
	GLfloat* t)//orig是向量起点，dir向量方向，v0,v1,v2为三角形三点坐标,如果向量交三角形，返回true,t的值为dir的因子     vec=orig+t*dir;其实t,u,v都可以不计算
{
	// E1
	GLfloat u, v;
	glm::vec3 E1 = v1 - v0;
	// E2
	glm::vec3 E2 = v2 - v0;
	// P
	//out_file << v1.x << " " << v1.y << " " << v1.z << endl;
	glm::vec3 P = glm::cross(dir, E2);//dir,E2也许写反了
	// determinant
	GLfloat det = glm::dot(E1, P);
	//out_file << E1.x << " " << E1.y << " " << E1.z <<" " << P.x<<" " << P.y<<" " << P.z<< " "  << endl;
	// keep det > 0, modify T accordingly
	glm::vec3 T;
	if (det >0)
		T = orig - v0;
	else{
		T = v0 - orig;
		det = -det;
	}
	// If determinant is near zero, ray lies in plane of triangle
	if (det < 0.0001f){
		//out_file << "det<0.0001   det ="<<det << endl;
		return false;
	}
	// Calculate u and make sure u <= 1
	u = glm::dot(T, P);
	if (u < 0.0f || u > det)
		return false;
	// Q
	glm::vec3 Q = glm::cross(T, E1);
	// Calculate v and make sure u + v <= 1
	v = glm::dot(dir, Q);
	if (v < 0.0f || u + v > det){
		//out_file << "v<0.0f" << endl;
		return false;
	}
	// Calculate t, scale parameters, ray intersects triangle
	*t = glm::dot(E2, Q);
	float fInvDet = 1.0f / det;
	*t *= fInvDet;
	u *= fInvDet;
	v *= fInvDet;
	return true;
}



bool interPointIndexLess(interPointIndex &a, interPointIndex &b){
	return a.distance < b.distance;
}

glm::vec3 computeNormal//计算三角形面的法向量
(
glm::vec3 const & a,
glm::vec3 const & b,
glm::vec3 const & c
)
{
	return glm::normalize(glm::cross(c - a, b - a));
}




float downx, downy;
int buttonState;
GLfloat Angle = 0.0, Angle2=0,Distance = 3.0;

//GLuint lightVAO, lightVBO, lightEBO, VBOs[NUMBEROFMODELS], VAOs[NUMBEROFMODELS], EBOs[NUMBEROFMODELS];
GLint numberOfPoints[NUMBEROFMODELS], numberOfTriangles[NUMBEROFMODELS], numberOfEdges[NUMBEROFMODELS];
GLint lnumberOfPoints, lnumberOfTriangles, lnumberOfEdges;
//GLfloat buffer1[100000], buffer2[100000];
//GLuint buffer3[100000], buffer4[100000];
//GLfloat *vertexs[NUMBEROFMODELS] = { buffer1, buffer2 };
//GLuint *indexs[NUMBEROFMODELS] = { buffer3, buffer4 };
GLfloat *vertexs[NUMBEROFMODELS], *lightVertexs;
Point *PointArray[NUMBEROFMODELS];
Index *IndexArray[NUMBEROFMODELS];
Triangle *Triangles[NUMBEROFMODELS];
GLuint *indexs[NUMBEROFMODELS], *lightIndexs;
char *fileNames[NUMBEROFMODELS] = {
	"boxTrace.off", "CupTrace.off","cup2Trace.off" };

glm::vec3 Colors[NUMBEROFMODELS] = { glm::vec3(1.0f, 1.0f, 0.0f),glm::vec3(0.0f,0.5f,0.5f),glm::vec3(1.0f,1.0f,1.0f)};


//====model
//glm::vec3 Positions[NUMBEROFMODELS] = { glm::vec3(0.0f, 0.0f, 0.0f) };
//glm::vec3  RotatePointers[NUMBEROFMODELS] = { glm::vec3(0.0f, 1.0f, 0.0f)};
//GLfloat RotateRadians[NUMBEROFMODELS] = { glm::radians(0.0f)};
//glm::vec3 Scales[NUMBEROFMODELS] = { glm::vec3(1.0f, 1.0f, 1.0f)};
//==========

//=====camera
bool keys[1024];
GLfloat lastX = 400, lastY = 300, yaw, pitch;
bool firstMouse = true;
bool LeftPress = false;
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 2.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
//=========

//====light
glm::vec3 lightPos(2.0f, 0.0f, 0.0f);
glm::vec3 lightSca(0.1f, 0.1f, 0.1f);
glm::vec3 lightRotPo(0.0f, 0.0f, 1.0f);
glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
GLfloat lightRoRadian = glm::radians(0.0f);
//=======

//======
GLfloat ambient=0.15f,specularStrength = 1.0f,reflectAbsolve=0.5f;



struct RGB{
	GLfloat R, G, B;
};

RGB Pixels[HEIGHT][WEIGHT];
/*


H



(0,0)        W

*/






class Ray3{//光线类
public:
	glm::vec3 origin, direction;
	vector<interPointIndex> inters;
	Ray3(glm::vec3 a, glm::vec3 b){
		origin = a;
		direction = b;
	}
	Ray3(Ray3 &a){
		origin = a.origin;
		direction = a.direction;
		inters = a.inters;
	}
	glm::vec3 getPoint(GLfloat t){//计算O+t*dir
		return origin + t*direction;
	}
	
};


class Camera{//暂时要求长宽相等
	glm::vec3 eye, front, refUp, right, up;
	GLfloat fovScale, fov;
public:
	Camera(glm::vec3 eye, glm::vec3 front, glm::vec3 up, GLfloat fov){//fov为视角（角度）
		this->eye = eye;
		this->front = front;
		this->refUp = up;
		this->fov = fov;
		initialize();
	}
	void initialize(){//当位置更新后需要调用这个
		this->right = glm::cross(front, refUp);
		this->up = glm::cross(right, front);
		this->fovScale = glm::tan(fov*0.5*PI / 180) * 2;//不知是否需要转化成角度
	}
	Ray3 generateRay(GLfloat x, GLfloat y){
		glm::vec3 r = right*((x - 0.5f)*fovScale);
		glm::vec3 u = up*((y - 0.5f)*fovScale);
		return Ray3(eye, glm::normalize(front + r + u));//返回对应光线
	}
};
Camera camera = Camera(cameraPos, cameraFront, cameraUp, 90);//500*500的窗口

class Sphere{//包围球，用于简单碰撞检测
	GLfloat sqrRadius;
public:
	GLfloat radius;
	glm::vec3 center;
	Sphere(glm::vec3 a, GLfloat b){ center = a; radius = b; sqrRadius = b*b; }
	bool isInter(Ray3 &ray){
		glm::vec3 v = ray.origin - center;
		GLfloat Ddotv = glm::dot(ray.direction, v);
		if (Ddotv <= 0){
			GLfloat a0 = glm::length(v)*glm::length(v) - sqrRadius;
			GLfloat discr = Ddotv*Ddotv - a0;
			if (discr >= 0){
				return true;
			}
		}
		return false;
	}
};
Sphere containers[NUMBEROFMODELS] = { Sphere(glm::vec3(0.0f, -1.5f, 0.0f), 1.0f) ,Sphere(glm::vec3(0.0f,1.0f,0.0f),1.0f),Sphere(glm::vec3(0.0f,0.0f,-1.0f),1.0f)};


interPointIndex getInterPointIndex(Ray3 &ray){//将光线和所有模型求交，并返回第一个相交面的下标，如果没有则下标为-1
	//out_file << ray.direction.x << " " << ray.direction.y << " " << ray.direction.z << " " << endl;
	if (containers[0].isInter(ray)){//如果与盒子相交才继续计算
		for (int i = 0; i < NUMBEROFMODELS; i++){
			for (int j = 0; j < numberOfTriangles[i]; j++){//遍历每个模型的所有面
				//int xindex = IndexArray[i][j].xindex, yindex = IndexArray[i][j].yindex, zindex = IndexArray[i][j].zindex;
				//cout << "i=" << i << " j=" << j << endl;
				GLfloat distance;
				if (IntersectTriangle(ray.origin, ray.direction,
					Triangles[i][j].P1, Triangles[i][j].P2, Triangles[i][j].P3, &distance)){//release下三角形坐标错了
					interPointIndex temp;
					temp.model = i;
					temp.triangle = j;
					temp.distance = distance*glm::length(ray.direction);
					//cout << temp.distance << endl;
					ray.inters.push_back(temp);//将相交面对应的模型和三角形下标以及距离保存起来。
				}
			}
		}
		interPointIndex temp;
		if (ray.inters.empty())
		{
			temp.model = -1;
			return temp;
		}
		else{
			temp = *min_element(ray.inters.begin(), ray.inters.end(), interPointIndexLess);
			return temp;
		}
	}
	interPointIndex temp;
	temp.model = -1;
	return temp;
}



glm::vec3 recurseColor(Ray3 &ray, int level){
	if (level == 3)
		return glm::vec3(0.0f);
	else{
		for (int i = 0; i < NUMBEROFMODELS; i++){
			if (containers[i].isInter(ray)){//如果与盒子相交才继续计算
				for (int j = 0; j < numberOfTriangles[i]; j++){//遍历每个模型的所有面
					//int xindex = IndexArray[i][j].xindex, yindex = IndexArray[i][j].yindex, zindex = IndexArray[i][j].zindex;
					//cout << "i=" << i << " j=" << j << endl;
					GLfloat distance;
					if (IntersectTriangle(ray.origin, ray.direction,
						Triangles[i][j].P1, Triangles[i][j].P2, Triangles[i][j].P3, &distance)){//release下三角形坐标错了
						interPointIndex temp;
						temp.model = i;
						temp.triangle = j;
						temp.distance = distance*glm::length(ray.direction);
						//cout << temp.distance << endl;
						ray.inters.push_back(temp);//将相交面对应的模型和三角形下标以及距离保存起来。
					}
				}
			}
		}
		if (ray.inters.empty())//如果没有相交的面，返回黑色
			return glm::vec3(0.0f);
		interPointIndex temp = *min_element(ray.inters.begin(), ray.inters.end(), interPointIndexLess);//取的相交面下标
		//==========反射光
		glm::vec3 viewDir = glm::normalize(ray.direction);
		glm::vec3 reflectDir = glm::reflect(viewDir, Triangles[temp.model][temp.triangle].norm);
		Ray3 reRay(Triangles[temp.model][temp.triangle].P1, reflectDir);
		glm::vec3 reflectRe = recurseColor(reRay, level + 1);//反射角度来的光
		//========计算漫反射和镜面反射
		glm::vec3 lightDir = glm::normalize(Triangles[temp.model][temp.triangle].P1 - lightPos);//任取平面一点
		interPointIndex toLight = getInterPointIndex(Ray3(Triangles[temp.model][temp.triangle].P1, lightDir));
		GLfloat diff,spec;
		if (toLight.model == -1){
			diff = glm::max(glm::dot(Triangles[temp.model][temp.triangle].norm, lightDir), 0.0f);//光照下漫反射因子
			spec = glm::pow(glm::max(glm::dot(reflectDir, -lightDir), 0.0f), 8);//lightDir可能反了
		}
		else{
			diff = glm::max(glm::dot(Triangles[temp.model][temp.triangle].norm, lightDir), 0.0f) / 3.0f;//被物体挡住的话，因子除3
			spec = 0.0f;
		}
		glm::vec3 result = Colors[temp.model] * (diff+ambient+spec) + reflectRe*reflectAbsolve;
		return result;
	}
}

void rayTrace(){
	cout << clock() << "ms" << endl;
	GLfloat delx = 1.0f / HEIGHT, dely = 1.0f / WEIGHT;
	for (int i = 0; i < HEIGHT; i++){
		//y = 0.0f;
		for (int j = 0; j < WEIGHT; j++){
			Ray3 ray = camera.generateRay((GLfloat)i/HEIGHT, (GLfloat)j/WEIGHT);
			//y += dely;
			/*interPointIndex temp = getInterPointIndex(ray);
			if (temp.model != -1)//光线碰到了物体,则要追踪光线取的颜色
			{
				Pixels[i][j].R = 1.0f*temp.distance / 10.0f;
				Pixels[i][j].G = 1.0f*temp.distance / 10.0f;
				Pixels[i][j].B = 1.0f*temp.distance / 10.0f;
			}
			else{
				Pixels[i][j].R = 0.0f;
				Pixels[i][j].G = 0.0f;
				Pixels[i][j].B = 0.0f;
			}*/
			glm::vec3 color = recurseColor(ray, 1);
			//out_file << color.x << " " << color.y << " " << color.z << endl;
			Pixels[i][j].R = color.x;
			Pixels[i][j].G = color.y;
			Pixels[i][j].B = color.z;
		}
		//x += delx;
	}
	cout << clock() << "ms" << endl;
}

void display(){
	//cout << Count << endl;
	//cout << Count2 << endl;
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0, 0.0, 0.0);
	glDrawPixels(WEIGHT, HEIGHT, GL_RGB, GL_FLOAT, (void *)Pixels);
	//glDrawPixels(500, 500, GL_RGB, GL_FLOAT, (void *)Pixels);
	//glEnd();
	glFlush();
}




void reshape(int w, int h){
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 20.0);
	//glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);//正投影
	glMatrixMode(GL_MODELVIEW);
}


void key_callback(int key, int x, int y)
{
	GLfloat cameraSpeed = 0.05f;
	if (key == GLUT_KEY_UP)
		cameraPos += cameraSpeed * cameraFront;
	if (key == GLUT_KEY_DOWN)
		cameraPos -= cameraSpeed * cameraFront;
	if (key == GLUT_KEY_LEFT)
		cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
	if (key == GLUT_KEY_RIGHT)
		cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
}


void init(){
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glShadeModel(GL_SMOOTH);
	//glEnable(GL_DEPTH_TEST);
	glLoadIdentity();
	gluLookAt(Distance*sin(Angle)*cos(Angle2), Distance*sin(Angle2), Distance*cos(Angle)*cos(Angle2), 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	//======
	
}

void readFile(char *fileName, int index){
	char buffer[20];
	fstream in_file(fileName, ios::in);
	if (!in_file){
		cout << "no such file" << endl;
		exit(-1);
	}
	char temp[10];
	in_file.getline(temp, 10);
	//cout << temp << endl;
	if (strcmp(temp, "OFF") != 0)
	{
		cout << "This file is not an OFF file" << endl;
		exit(-1);
	}
	in_file >> numberOfPoints[index] >> numberOfTriangles[index] >> numberOfEdges[index];
	//cout << numberOfPoints << "  " << numberOfTriangles << endl;
	vertexs[index] = new GLfloat[numberOfPoints[index] * 6];
	for (int i = 0; i < numberOfPoints[index] * 6; i++)//初始化
		vertexs[index][i] = 0.0f;
	PointArray[index] = (Point *)(vertexs[index]);
	//indexs[index] = new GLuint[numberOfTriangles[index] * 3];
	Triangles[index] = new Triangle[numberOfTriangles[index]];
	//IndexArray[index] = (Index *)(indexs[index]);
	/*for (int i = 0; i < 6 * numberOfPoints[index]; i++){
		in_file >> vertexs[index][i++] >> vertexs[index][i++] >> vertexs[index][i++];
		i += 2;
		//cout<< vdata[i][0] << "   " << vdata[i][1] << "   " << vdata[i][2] << endl;
	}*/
	for (int i = 0; i < numberOfPoints[index]; i++){//release版本下用上面的会发生读取错误？导致坐标的y,z没有读进来
		in_file >> PointArray[index][i].x >> PointArray[index][i].y >> PointArray[index][i].z;
	}
	int count = 0;
	for (int i = 0; i < numberOfTriangles[index]; i++){//还要计算法向量
		int temp;
		in_file >> temp;//丢弃第一个数，因为所有面都是三角形
		//in_file >> indexs[index][count++] >> indexs[index][count++] >> indexs[index][count++];
		GLuint a, b, c;
		in_file >> a >> b >> c;
		glm::vec3 P1 = glm::vec3(PointArray[index][a].x, PointArray[index][a].y, PointArray[index][a].z);
		glm::vec3 P2 = glm::vec3(PointArray[index][b].x, PointArray[index][b].y, PointArray[index][b].z);
		glm::vec3 P3 = glm::vec3(PointArray[index][c].x, PointArray[index][c].y, PointArray[index][c].z);//三角形三个点坐标
		//out_file << PointArray[index][a].x << " " << PointArray[index][a].y << " " << PointArray[index][a].z << endl;//PointArray出错
		glm::vec3 norm = computeNormal(P1, P2, P3);
		Triangles[index][i].P1 = P1;
		Triangles[index][i].P2 = P2;
		Triangles[index][i].P3 = P3;
		Triangles[index][i].norm = norm;
	}
}

// Determine whether a ray intersect with a triangle
// Parameters
// orig: origin of the ray
// dir: direction of the ray
// v0, v1, v2: vertices of triangle
// t(out): weight of the intersection for the ray
// u(out), v(out): barycentric coordinate of intersection




int main(int argc, char ** argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize(WEIGHT, HEIGHT);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);
	init();
	for (int i = 0; i < NUMBEROFMODELS; i++)
		readFile(fileNames[i], i);
	//computeNorm();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	rayTrace();
	glutMainLoop();
	return 0;
}