/*
 * util.cpp
 *
 *	Implementation of functions in util.h
 *
 *  Created on: Mar 9, 2010
 *      Author: drogers
 */

#include "util.h"

#include <cstdlib>
#include <cmath>

using namespace std;

///******* Vector stuff

/**
 * Returns magnitude or length of vector;
 */
GLfloat magnitude(GLfloat vec[3]) {
	return (GLfloat) sqrt( pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
}

/**
 * Normalizes vector, ie turns it into a unit vector with the same direction.
 */
void normalize(GLfloat vec[3]) {
	GLfloat mag = magnitude(vec);
	if(mag == 0) {
		vec[0] = vec[1] = vec[2] = 0;
	} else {
		for(int i=0; i<3; i++) {
			vec[i] /= mag;
		}
	}
}

/**
 * Puts the crossproduct of left X right into out.
 */
void cross(GLfloat out[3], GLfloat left[3], GLfloat right[3]) {
	out[0] = left[1]*right[2] - left[2]*right[1]; // y*v.z - z*v.y
	out[1] = left[2]*right[0] - left[0]*right[2];// z*v.x - x*v.z
	out[2] = left[0]*right[1] - left[1]*right[0];// x*v.y - y*v.x
}

/**
 * Copies values from copy_from to copy_to.
 */
void copypv(GLfloat copy_to[3], GLfloat copy_from[3]) {
	for (int i = 0; i < 3; ++i) {
		copy_to[i] = copy_from[i];
	}
}

/**
 * Set values in vec to zero.
 */
void zero(GLfloat vec[3]) {
	vec[0] = 0;
	vec[1] = 0;
	vec[2] = 0;
}

/**
 * Dot product for vectors, vec_sz defaults to 3, specify for other dim.
 */
GLfloat dot(GLfloat *v1, GLfloat *v2, int vec_sz) {
	GLfloat dot = 0;
	for (int i = 0; i < vec_sz; ++i) {
		dot += v1[i]*v2[i];
	}
	return dot;
}

////********************* Matrix Stuff ***********************////

void set_to_ident(GLfloat m[16]) {
	for (int i = 0; i < 16; ++i) {
		m[i] = 0;
	}
	m[0] = 1;
	m[5] = 1;
	m[10] = 1;
	m[15] = 1;
}

/**
 * Copies values from copy_from to copy_to.
 */
void copym(GLfloat copy_to[16], GLfloat copy_from[16]) {
	for (int i = 0; i < 16; ++i) {
		copy_to[i] = copy_from[i];
	}
}

/**
 * Multiply left * right and put product into out.
 * ie: out <-- left right
 */
void mult_matrixf(GLfloat out[16], GLfloat left[16], GLfloat right[16]) {
	GLfloat rows[4][4], cols[4][4];
	for (int i = 0; i < 4; ++i) {

		for (int j = 0; j < 4; ++j) {
			rows[i][j] = left[i*4 + j];
			cols[i][j] = right[i + 4*j];
		}
	}

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			out[i*4 + j] = dot(rows[i], cols[j], 4);
		}
	}

}
