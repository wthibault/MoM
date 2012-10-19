/*
 * util.h
 *
 *	Utility stuff with no object oriented stuff.
 *
 *  Created on: Mar 9, 2010
 *      Author: drogers
 */

#ifndef UTIL_H_
#define UTIL_H_


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

const GLfloat PI = 3.14159265;

///******* Vector stuff

/**
 * Returns magnitude or length of vector;
 */
GLfloat magnitude(GLfloat vec[3]);

/**
 * Normalizes vector, ie turns it into a unit vector with the same direction.
 */
void normalize(GLfloat vec[3]);

/**
 * Puts the crossproduct of left X right into out.
 */
void cross(GLfloat out[3], GLfloat left[3], GLfloat right[3]);

/**
 * Copies values from copy_from to copy_to.
 */
void copypv(GLfloat copy_to[3], GLfloat copy_from[3]);

/**
 * Dot product for vectors, vec_sz defaults to 3, specify for other dim.
 */
GLfloat dot(GLfloat *v1, GLfloat *v2, int vec_sz=3);

////************************* Matrix Stuff **********************////
/**
 * Sets column-order (ie opengl) matrix m to identity matrix.
 */
void set_to_ident(GLfloat m[16]);

/**
 * Copies values from copy_from to copy_to.
 */
void copym(GLfloat copy_to[16], GLfloat copy_from[16]);

/**
 * Multiply left * right and put product into out.
 * ie: out <-- left right
 * Does row major order.
 */
void mult_matrixf(GLfloat out[16], GLfloat left[16], GLfloat right[16]);



#endif /* UTIL_H_ */
