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
