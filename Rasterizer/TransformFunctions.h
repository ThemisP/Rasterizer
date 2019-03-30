#pragma once
#include <glm/glm.hpp>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::mat4;

class TransformFunctions {
public:
	static void RotateAroundX(mat4 &matrix, float angle);
	static void RotateAroundY(mat4 &matrix, float angle);
private:
	static void Rotate(mat3 rotation, mat4 &matrix);
};

