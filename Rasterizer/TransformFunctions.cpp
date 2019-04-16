#include "TransformFunctions.h"

void TransformFunctions::RotateAroundX(mat4 &matrix, float angle) {
	vec3 x0(1, 0, 0);
	vec3 x1(0, cosf(angle), sinf(angle));
	vec3 x2(0, -sinf(angle), cosf(angle));

	Rotate(mat3(x0, x1, x2), matrix);
}

void TransformFunctions::RotateAroundY(mat4 &matrix, float angle) {
	vec3 x0(cosf(angle), 0, -sinf(angle));
	vec3 x1(0, 1, 0);
	vec3 x2(sinf(angle), 0, cosf(angle));

	Rotate(mat3(x0, x1, x2), matrix);
}

mat4 TransformFunctions::LookAt(vec3 cameraPos, vec3 target, vec3 up) {
	vec3 cameraDir = glm::normalize(cameraPos - target);
	vec3 rightDir = glm::normalize(glm::cross(up, cameraDir));
	vec3 upDir = glm::normalize(glm::cross(cameraDir, rightDir));

	mat4 matrix1(rightDir[0], rightDir[1], rightDir[2], 0,
				 upDir[0], upDir[1], upDir[2], 0,
				 cameraDir[0], cameraDir[1], cameraDir[2], 0,
				 0, 0, 0, 1);
	mat4 matrix2(1, 0, 0, -cameraPos[0],
				 0, 1, 0, -cameraPos[1],
				 0, 0, 1, -cameraPos[2],
				 0, 0, 0, 1);
	return matrix1 * matrix2;
}

void TransformFunctions::Rotate(mat3 rotation, mat4 &matrix) {
	mat3 rotationExtract(matrix[0][0], matrix[0][1], matrix[0][2],
		matrix[1][0], matrix[1][1], matrix[1][2],
		matrix[2][0], matrix[2][1], matrix[2][2]);

	rotationExtract = rotation * rotationExtract;

	matrix[0][0] = rotationExtract[0][0];
	matrix[0][1] = rotationExtract[0][1];
	matrix[0][2] = rotationExtract[0][2];
	matrix[1][0] = rotationExtract[1][0];
	matrix[1][1] = rotationExtract[1][1];
	matrix[1][2] = rotationExtract[1][2];
	matrix[2][0] = rotationExtract[2][0];
	matrix[2][1] = rotationExtract[2][1];
	matrix[2][2] = rotationExtract[2][2];
}
