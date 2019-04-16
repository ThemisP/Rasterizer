#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include "TransformFunctions.h"
#include <stdint.h>
#include <limits.h>

using namespace std;
using glm::ivec2;
using glm::vec2;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

//#define SCREEN_WIDTH 320
//#define SCREEN_HEIGHT 256
#define SCREEN_WIDTH 1024
#define SCREEN_HEIGHT 720
#define FULLSCREEN_MODE false

struct Pixel
{
	int x;
	int y;
	float zinv;
	vec4 pos3d;
};

struct Vertex {
	vec4 pos;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
vector<Triangle> testScene;
vec4 cameraPos(0, 0, -3.001, 1);
mat4 cameraTransform;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH] = { 0 };
vec4 lightPos(0, -0.5, -0.7, 1);
vec3 lightPower = 14.0f*vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);
mat4 view;
float time = 0.0f;

vec4 currentNormal;
vec3 currentReflectance;

bool Update();
void Draw(screen* screen);
void VertexShader(const Vertex& v, Pixel& p);
void PixelShader(screen* screen, const Pixel& p, vec3 color);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
//void DrawLineSDL(screen* screen, Pixel a, Pixel b, vec3 color);
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec3 color);
void DrawPolygonEdges(screen* screen, const vector<Vertex>& vertices);
void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 color);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);

vector<Vertex> ClipTriangle(vector<Vertex> vertices);
vector<Vertex> ClipTop(vector<Vertex> vertices);
vector<Vertex> ClipBot(vector<Vertex> vertices);
vector<Vertex> ClipRight(vector<Vertex> vertices);
vector<Vertex> ClipLeft(vector<Vertex> vertices);
vector<Vertex> ClipFront(vector<Vertex> vertices);
vector<Vertex> ClipBack(vector<Vertex> vertices);


int main(int argc, char* argv[])
{
	screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
	LoadTestModel(testScene);

	vec4 matr1(1, 0, 0, 0);
	vec4 matr2(0, 1, 0, 0);
	vec4 matr3(0, 0, 1, 0);
	vec4 matr4(0, 0, 0, 1);
	cameraTransform = mat4(matr1, matr2, matr3, matr4);

	while (Update())
	{
		Draw(screen);
		SDL_Renderframe(screen);
	}

	SDL_SaveImage(screen, "screenshot.bmp");

	KillSDL(screen);
	return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
	/* Clear buffer */
	memset(screen->buffer, 0, screen->height*screen->width * sizeof(uint32_t));
	memset(depthBuffer, 0, SCREEN_WIDTH*SCREEN_HEIGHT * sizeof(float));

	for (uint32_t i = 0; i < testScene.size(); ++i)
	{
		vector<Vertex> vertices(3);
		vec4 v0 = testScene[i].v0;
		vec4 v1 = testScene[i].v1;
		vec4 v2 = testScene[i].v2;
		vertices[0].pos = v0;
		vertices[1].pos = v1;
		vertices[2].pos = v2;

		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 norm = glm::cross(e1, e2);
		vec4 normal = vec4(glm::normalize(norm), 1);
		currentNormal = normal;

		vector<Vertex> clippingVertices;
		clippingVertices = ClipTriangle(vertices);
		int size = clippingVertices.size();
		if (size > 3) {
			for (int i = 2; i < size; i++) {
				vector<Vertex> vertClip;
				vertClip.push_back(clippingVertices[i-1]);
				vertClip.push_back(clippingVertices[i]);
				vertClip.push_back(clippingVertices[(i + 1) % size]);
				DrawPolygon(screen, vertClip, testScene[i].color);
			}
		} else {
			DrawPolygon(screen, clippingVertices, testScene[i].color);
		}

		
		//DrawPolygonEdges(screen, vertices);
	}

	//for (int i = 0; i < SCREEN_HEIGHT; i++) {
	//	for (int j = 0; j < SCREEN_WIDTH; j++) {
	//		float value = depthBuffer[i][j];
	//		vec3 color(value, value, value);
	//		color *= 0.5;
	//		//PutPixelSDL(screen, j, i, color);
	//	}
	//}

}

/*Place updates of parameters here*/
bool Update()
{
	static int t = SDL_GetTicks();
	/* Compute frame time */
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;

	SDL_Event e;
	while (SDL_PollEvent(&e))
	{
		
		if (e.type == SDL_QUIT)
		{
			return false;
		} else
			if (e.type == SDL_KEYDOWN)
			{
				/*float radius = 1.f;
				float camX = sin(time) * radius;
				float camZ = cos(time) * radius;
				time += 0.05f;
				glm::mat4 view;
				view = TransformFunctions::LookAt(glm::vec3(camX, 0.0, -camZ), glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 1.0, 0.0));
				cameraPos = cameraPos * view;*/
				int key_code = e.key.keysym.sym;
				switch (key_code)
				{
				case SDLK_UP:
					/* Move camera forward */
					TransformFunctions::RotateAroundX(cameraTransform, 0.1f);
					break;
				case SDLK_DOWN:
					/* Move camera backwards */
					TransformFunctions::RotateAroundX(cameraTransform, -0.1f);
					break;
				case SDLK_LEFT:
					/* Move camera left */
					TransformFunctions::RotateAroundY(cameraTransform, -0.1f);
					break;
				case SDLK_RIGHT:
					/* Move camera right */
					TransformFunctions::RotateAroundY(cameraTransform, 0.1f);
					break;
				case SDLK_w:
					lightPos.z += 0.5;
					break;
				case SDLK_s:
					lightPos.z -= 0.5;
					break;
				case SDLK_a:
					lightPos.x -= 0.5;
					break;
				case SDLK_d:
					lightPos.x += 0.5;
					break;
				case SDLK_i:
					cameraPos.z += 0.2;
					break;
				case SDLK_k:
					cameraPos.z -= 0.2;
					break;
				case SDLK_j:
					cameraPos.x -= 0.2;
					break;
				case SDLK_l:
					cameraPos.x += 0.2;
					break;
				case SDLK_ESCAPE:
					/* Move camera quit */
					return false;
				}
			}
	}
	return true;
}

void DrawPolygonEdges(screen* screen, const vector<Vertex>& vertices)
{
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<Pixel> projectedVertices(V);
	for (int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}
	// Loop over all vertices and draw the edge from it to the next vertex:
	for (int i = 0; i < V; ++i)
	{
		int j = (i + 1) % V; // The next vertex
		vec3 color(1, 1, 1);
		//DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	}
}

void VertexShader(const Vertex& v, Pixel& p) {
	vec4 point = v.pos -cameraPos * cameraTransform;
	p.x = round(SCREEN_HEIGHT * (point.x / point.z) + SCREEN_WIDTH * 0.5);
	p.y = round(SCREEN_HEIGHT * (point.y / point.z) + SCREEN_HEIGHT * 0.5);
	p.zinv = (float)(1 / point.z);
	p.pos3d = v.pos;


	/*vec3 rVec = v.pos - lightPos;
	vec3 D = (lightPower*(fmaxf(glm::dot(rVec, (vec3)v.normal), 0)));
	float denominator = 4 * M_PI*(glm::dot(rVec,rVec));
	D.x = D.x / denominator;
	D.y = D.y / denominator;
	D.z = D.z / denominator;
	vec3 reflectance = v.reflectance*(D + indirectLightPowerPerArea);
	p.illumination = reflectance;*/
}

void PixelShader(screen* screen, const Pixel& p, vec3 color) {
	int x = p.x;
	int y = p.y;
	if (p.zinv > depthBuffer[y][x])
	{
		vec3 rVec = p.pos3d - lightPos;
		vec3 D = (lightPower*(fmaxf(glm::dot(rVec, (vec3)currentNormal), 0)));
		float denominator = 4 * M_PI*(glm::dot(rVec, rVec));
		D.x = D.x / denominator;
		D.y = D.y / denominator;
		D.z = D.z / denominator;
		vec3 reflectance = color*(D + indirectLightPowerPerArea);
		depthBuffer[y][x] = p.zinv;
		PutPixelSDL(screen, x, y, reflectance);
	}
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result) {
	int N = result.size();
	vec2 step = vec2(b - a) / float(std::fmax(N - 1, 1));
	vec2 current(a);
	for (int i = 0; i < N; ++i)
	{
		result[i] = round(current);
		current += step;
	}
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result) {
	int N = result.size();
	float x = (b.x - a.x) / float(std::fmaxf(N - 1, 1));
	float y = (b.y - a.y) / float(std::fmaxf(N - 1, 1));
	float d = sqrt((b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y));
	float m = sqrt(x*x + y*y);
	float stepInv = m/d;

	Pixel current(a);
	
	for (int i = 0; i < N; ++i)	{
		result[i] = current;
		current.x = round(a.x+i*x);
		current.y = round(a.y+i*y);
		current.zinv = (a.zinv*(1 - (i*stepInv)) + b.zinv*(i*stepInv));
		vec4 pos1 = a.pos3d*a.zinv*(1 - (i*stepInv));
		vec4 pos2 = b.pos3d*b.zinv*(i*stepInv);
		current.pos3d = (float)(1 / current.zinv)*(pos1 + pos2);
		/*vec3 col1 = a.illumination*a.zinv*(1 - (i*stepInv));
		vec3 col2 = b.illumination*b.zinv*(i*stepInv);
		current.illumination = (1/current.zinv)*(col1+col2);*/
	}
}

void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec3 color) {
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for (int i = 0; i < V; ++i)
		VertexShader(vertices[i], vertexPixels[i]);
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(screen, leftPixels, rightPixels, color);
}

//
//void DrawLineSDL(screen* screen, Pixel a, Pixel b, vec3 color) {
//	ivec2 delta = glm::abs(a - b);
//	int pixels = glm::max(delta.x, delta.y) + 1;
//	vector<ivec2> results(pixels);
//	Interpolate(a, b, results);
//
//	for (int i = 0; i < results.size(); i++) {
//		PutPixelSDL(screen, results[i].x, results[i].y, color);
//	}
//}

void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 color) {
	for (int i = 0; i < leftPixels.size(); i++) {
		int length = rightPixels[i].x - leftPixels[i].x+1;

		vector<Pixel> results(length);
		Interpolate(leftPixels[i], rightPixels[i], results);
		for (Pixel pixel : results) {
			if ((pixel.x >= 0) && (pixel.x < SCREEN_WIDTH) && (pixel.y >= 0) && (pixel.y < SCREEN_HEIGHT)) {
				PixelShader(screen, pixel, color);
				/*if (depthBuffer[pixel.y][pixel.x] < pixel.zinv) {

					PutPixelSDL(screen, pixel.x, pixel.y, color);
					depthBuffer[pixel.y][pixel.x] = pixel.zinv;
				} */
				//PutPixelSDL(screen, pixel.x, pixel.y, color);
			}
		}
		
	}
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int ymin = + numeric_limits<int>::max();
	int ymax = - numeric_limits<int>::max();
	for (Pixel vertex : vertexPixels) {
		if (vertex.y < ymin)
			ymin = vertex.y;
		if (vertex.y > ymax)
			ymax = vertex.y;
	}

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	int Rows = ymax - ymin + 1;
	leftPixels.resize(Rows);
	rightPixels.resize(Rows);	

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < Rows; i++) {
		leftPixels[i].x = +numeric_limits<int>::max();
		leftPixels[i].y = ymin + i;
		rightPixels[i].x = -numeric_limits<int>::max();
		rightPixels[i].y = ymin + i;
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (int i = 0; i < vertexPixels.size(); ++i)
	{
		int j = (i + 1) % vertexPixels.size(); // The next vertex
		vec3 color(1, 1, 1);
		//DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
		float deltax = glm::abs(vertexPixels[i].x - vertexPixels[j].x);
		float deltay = glm::abs(vertexPixels[i].y - vertexPixels[j].y);
		int pixels = glm::max(deltax, deltay) + 1;
		vector<Pixel> results(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], results);
		for (Pixel result : results) {
			int loc = result.y - ymin;
			if (result.x < leftPixels[loc].x) {
				leftPixels[loc].x = result.x;
				leftPixels[loc].zinv = result.zinv;
				leftPixels[loc].pos3d = result.pos3d;
			}
			if (result.x > rightPixels[loc].x) {
				rightPixels[loc].x = result.x;
				rightPixels[loc].zinv = result.zinv;
				rightPixels[loc].pos3d = result.pos3d;
			}
		}
	}
}

vector<Vertex> ClipTriangle(vector<Vertex> vertices) {
	vector<Vertex> clipped = vertices;
	//clipped = ClipTop(clipped);
	//clipped = ClipBot(clipped);
	//clipped = ClipRight(clipped);
	//clipped = ClipLeft(clipped);
	//clipped = ClipFront(clipped);
	//clipped = ClipBack(clipped);
	return clipped;
}

vector<Vertex> ClipBack(vector<Vertex> vertices) {
	vector<Vertex> clipped;
	return clipped;
}

vector<Vertex> ClipFront(vector<Vertex> vertices) {
	vector<Vertex> clipped;
	return clipped;
}

vector<Vertex> ClipLeft(vector<Vertex> vertices) {
	vector<Vertex> clipped;
	return clipped;
}

vector<Vertex> ClipRight(vector<Vertex> vertices) {
	vector<Vertex> clipped;
	return clipped;
}

vector<Vertex> ClipBot(vector<Vertex> vertices) {
	vector<Vertex> clipped;
	float yMin = -SCREEN_HEIGHT * 0.5 - 100;
	int size = vertices.size();
	for (int i = 0; i < size; i++) {
		Vertex v1 = vertices[i];
		vec4 pos1 = v1.pos - cameraPos * cameraTransform;
		Vertex v2 = vertices[(i + 1) % size];
		vec4 pos2 = v2.pos - cameraPos * cameraTransform;
		float v1class = pos1.y - yMin;
		float v2class = pos2.y - yMin;

		if (v1class > 0 && v2class > 0) {
			clipped.push_back(v2);
		}

		if (v1class > 0 && v2class < 0) {
			float step = v1class / (v1class - v2class);
			Vertex v;
			v.pos = v1.pos + step * (v1.pos - v2.pos);
			clipped.push_back(v);
		}
		if (v2class > 0 && v1class < 0) {
			float step = v2class / (v2class - v1class);
			Vertex v;
			v.pos = v2.pos + step * (v2.pos - v1.pos);
			clipped.push_back(v);
			clipped.push_back(v2);
		}
	}
	return clipped;
}

vector<Vertex> ClipTop(vector<Vertex> vertices) {
	vector<Vertex> clipped;
	float yMax = SCREEN_HEIGHT + 100;
	int size = vertices.size();
	for (int i = 0; i < size; i++) {
		Vertex v1 = vertices[i];
		float y1 = round(SCREEN_HEIGHT * (v1.pos.y / v1.pos.z) + SCREEN_WIDTH * 0.5);
		Vertex v2 = vertices[(i + 1) % size];
		float y2 = round(SCREEN_HEIGHT * (v2.pos.y / v2.pos.z) + SCREEN_WIDTH * 0.5);
		float v1class = yMax - y1;
		float v2class = yMax - y2;

		if (v1class > 0 && v2class > 0) {
			clipped.push_back(v2);
		}

		if (v1class > 0 && v2class < 0) {
			float step = y1 / (y1 - y2);
			Vertex v;
			v.pos = v1.pos + step * (v1.pos - v2.pos);
			clipped.push_back(v);
		}
		if (v2class > 0 && v1class < 0) {
			float step = y2 / (y2 - y1);
			Vertex v;
			v.pos = v2.pos + step * (v2.pos - v1.pos);
			clipped.push_back(v);
			clipped.push_back(v2);
		}
	}
	return clipped;
}