/*
Vector Matrix Library

MIT License

Copyright (c) 2024 Dice

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef VML_H
#define VML_H

#define VML_VERSION_MAJOR 1
#define VML_VERSION_MINOR 0
#define VML_VERSION_PATCH 0

// SCALAR FUNCTIONS

static inline float vml_min(
	float a,
	float b
);

static inline float vml_max(
	float a,
	float b
);

static inline float vml_abs(
	float a
);

static inline float vml_sign(
	float a
);

static inline float vml_lerp(
	float a,
	float b,
	float t
);

// VEC2 FUNCTIONS

static inline void vml_vec2_neg(
	float a[2],
	float b[2]
);

static inline void vml_vec2_abs(
	float a[2],
	float b[2]
);

static inline void vml_vec2_mov(
	float a[2],
	float b[2]
);

static inline void vml_vec2_add(
	float a[2],
	float b[2],
	float c[2]
);

static inline void vml_vec2_sub(
	float a[2],
	float b[2],
	float c[2]
);

static inline void vml_vec2_mul(
	float a[2],
	float b[2],
	float c[2]
);

static inline void vml_vec2_div(
	float a[2],
	float b[2],
	float c[2]
);

static inline void vml_vec2_dot(
	float  a[2],
	float  b[2],
	float *c
);

static inline void vml_vec2_unit(
	float a[2],
	float b[2]
);

static inline void vml_vec2_mag(
	float  a[2],
	float *b
);

static inline void vml_vec2_lerp(
	float a[2],
	float b[2],
	float t,
	float c[2]
);

// VEC3 FUNCTIONS

static inline void vml_vec3_neg(
	float a[3],
	float b[3]
);

static inline void vml_vec3_abs(
	float a[3],
	float b[3]
);

static inline void vml_vec3_mov(
	float a[3],
	float b[3]
);

static inline void vml_vec3_add(
	float a[3],
	float b[3],
	float c[3]
);

static inline void vml_vec3_sub(
	float a[3],
	float b[3],
	float c[3]
);

static inline void vml_vec3_mul(
	float a[3],
	float b[3],
	float c[3]
);

static inline void vml_vec3_div(
	float a[3],
	float b[3],
	float c[3]
);

static inline void vml_vec3_dot(
	float  a[3],
	float  b[3],
	float *c
);

static inline void vml_vec3_cross(
	float a[3],
	float b[3],
	float c[3]
);

static inline void vml_vec3_unit(
	float a[3],
	float b[3]
);

static inline void vml_vec3_mag(
	float  a[3],
	float *b
);

static inline void vml_vec3_lerp(
	float a[3],
	float b[3],
	float t,
	float c[3]
);

// QUATERNION FUNCTIONS

void vml_quat_lerp(
	float a[4],
	float b[4],
	float t,
	float c[4]
);

// MAT3 FUNCTIONS

static inline void vml_mat3_id(
	float a[9]
);

static inline void vml_mat3_mov(
	float a[9],
	float b[9]
);

void vml_mat3_transpose(
	float a[9],
	float b[9]
);

void vml_mat3_vec3_mul(
	float a[9],
	float b[3],
	float c[3]
);

void vml_mat3_mul(
	float a[9],
	float b[9],
	float c[9]
);

void vml_mat3_inv(
	float a[9],
	float b[9]
);

void vml_mat3_scale_set(
	float a[9],
	float b[3],
	float c[9]
);

// MAT4 FUNCTIONS

static inline void vml_mat4_id(
	float a[16]
);

static inline void vml_mat4_mov(
	float a[16],
	float b[16]
);

void vml_mat4_transpose(
	float a[16],
	float b[16]
);

void vml_mat4_vec3_mul(
	float a[16],
	float b[3],
	float c[3]
);

void vml_mat4_mul(
	float a[16],
	float b[16],
	float c[16]
);

void vml_mat4_inv(
	float a[16],
	float b[16]
);

void vml_mat4_mat3_get(
	float a[16],
	float b[9]
);

void vml_mat4_euler_set(
	float a[16],
	float b[3],
	float c[16]
);

void vml_mat4_euler_get(
	float a[16],
	float b[3]
);

void vml_mat4_quat_set(
	float a[16],
	float b[4],
	float c[16]
);

void vml_mat4_quat_get(
	float a[16],
	float b[4]
);

void vml_mat4_axis_set(
	float a[16],
	float b[4],
	float c[16]
);

void vml_mat4_axis_get(
	float a[16],
	float b[4]
);

void vml_mat4_look_set(
	float a[16],
	float b[3],
	float c[3],
	float d[16]
);

void vml_mat4_look_get(
	float a[16],
	float b[3]
);

void vml_mat4_position_set(
	float a[16],
	float b[3],
	float c[16]
);

void vml_mat4_position_get(
	float a[16],
	float b[3]
);

void vml_mat4_scale_set(
	float a[16],
	float b[3],
	float c[16]
);

void vml_mat4_scale_get(
	float a[16],
	float b[3]
);

void vml_mat4_lerp(
	float a[16],
	float b[16],
	float t,
	float c[16]
);

void vml_mat4_perspective(
	float v, // Field of view
	float r, // Aspect ratio
	float n, // Near plane
	float f, // Far plane
	float a[16]
);

void vml_mat4_orthographic(
	float l, // Left
	float r, // Right
	float t, // Top
	float b, // Bottom
	float n, // Near
	float f, // Far
	float a[16]
);

#endif

/************************[IMPLEMENTATION BEGINS HERE]*************************/

#ifdef VML_IMPLEMENTATION
#ifndef VML_C
#define VML_C

#include <math.h>

// SCALAR FUNCTIONS

static inline float vml_min(
	float a,
	float b
) {
	return (a<b)?a:b;
}

static inline float vml_max(
	float a,
	float b
) {
	return (a>b)?a:b;
}

static inline float vml_abs(
	float a
) {
	return (a<0)?-a:a;
}

static inline float vml_sign(
	float a
) {
	if (a==0) return 0;

	return (a<0)?-1:1;
}

static inline float vml_lerp(
	float a,
	float b,
	float t
) {
	return a*(1-t)+b*t;
}

// VEC2 FUNCTIONS

static inline void vml_vec2_neg(
	float a[2],
	float b[2]
) {
	b[0] = -a[0];
	b[1] = -a[1];
}

static inline void vml_vec2_abs(
	float a[2],
	float b[2]
) {
	b[0] = vml_abs(a[0]);
	b[1] = vml_abs(a[1]);
}

static inline void vml_vec2_mov(
	float a[2],
	float b[2]
) {
	b[0] = a[0];
	b[1] = a[1];
}

static inline void vml_vec2_add(
	float a[2],
	float b[2],
	float c[2]
) {
	c[0] = a[0]+b[0];
	c[1] = a[1]+b[1];
}

static inline void vml_vec2_sub(
	float a[2],
	float b[2],
	float c[2]
) {
	c[0] = a[0]-b[0];
	c[1] = a[1]-b[1];
}

static inline void vml_vec2_mul(
	float a[2],
	float b[2],
	float c[2]
) {
	c[0] = a[0]*b[0];
	c[1] = a[1]*b[1];
}

static inline void vml_vec2_div(
	float a[2],
	float b[2],
	float c[2]
) {
	c[0] = a[0]/b[0];
	c[1] = a[1]/b[1];
}

static inline void vml_vec2_dot(
	float  a[2],
	float  b[2],
	float *c
) {
	*c = (a[0]*b[0])+(a[1]*b[1]);
}

static inline void vml_vec2_unit(
	float a[2],
	float b[2]
) {
	float m=sqrtf(
		a[0]*a[0]+
		a[1]*a[1]
	);

	if (m==0) {
		b[0] = 0;
		b[1] = 0;

		return;
	}

	b[0] = a[0]/m;
	b[1] = a[1]/m;
}

static inline void vml_vec2_mag(
	float  a[2],
	float *b
) {
	*b=sqrtf(
		a[0]*a[0]+
		a[1]*a[1]
	);
}

static inline void vml_vec2_lerp(
	float a[2],
	float b[2],
	float t,
	float c[2]
) {
	float it = 1-t;

	c[0] = vml_lerp(a[0],b[0],t);
	c[1] = vml_lerp(a[1],b[1],t);
}

// VEC3 FUNCTIONS

static inline void vml_vec3_neg(
	float a[3],
	float b[3]
) {
	b[0] = -a[0];
	b[1] = -a[1];
	b[2] = -a[2];
}

static inline void vml_vec3_abs(
	float a[3],
	float b[3]
) {
	b[0] = vml_abs(a[0]);
	b[1] = vml_abs(a[1]);
	b[2] = vml_abs(a[2]);
}

static inline void vml_vec3_mov(
	float a[3],
	float b[3]
) {
	b[0] = a[0];
	b[1] = a[1];
	b[2] = a[2];
}

static inline void vml_vec3_add(
	float a[3],
	float b[3],
	float c[3]
) {
	c[0] = a[0]+b[0];
	c[1] = a[1]+b[1];
	c[2] = a[2]+b[2];
}

static inline void vml_vec3_sub(
	float a[3],
	float b[3],
	float c[3]
) {
	c[0] = a[0]-b[0];
	c[1] = a[1]-b[1];
	c[2] = a[2]-b[2];
}

static inline void vml_vec3_mul(
	float a[3],
	float b[3],
	float c[3]
) {
	c[0] = a[0]*b[0];
	c[1] = a[1]*b[1];
	c[2] = a[2]*b[2];
}

static inline void vml_vec3_div(
	float a[3],
	float b[3],
	float c[3]
) {
	c[0] = a[0]/b[0];
	c[1] = a[1]/b[1];
	c[2] = a[2]/b[2];
}

static inline void vml_vec3_dot(
	float  a[3],
	float  b[3],
	float *c
) {
	*c = (a[0]*b[0])+(a[1]*b[1])+(a[2]*b[2]);
}

static inline void vml_vec3_cross(
	float a[3],
	float b[3],
	float c[3]
) {
	float r0 = a[1]*b[2]-a[2]*b[1];
	float r1 = a[2]*b[0]-a[0]*b[2];
	float r2 = a[0]*b[1]-a[1]*b[0];

	c[0] = r0;
	c[1] = r1;
	c[2] = r2;
}

static inline void vml_vec3_unit(
	float a[3],
	float b[3]
) {
	float m=sqrtf(
		a[0]*a[0]+
		a[1]*a[1]+
		a[2]*a[2]
	);

	if (m==0) {
		b[0] = 0;
		b[1] = 0;
		b[2] = 0;

		return;
	}

	b[0] = a[0]/m;
	b[1] = a[1]/m;
	b[2] = a[2]/m;
}

static inline void vml_vec3_mag(
	float  a[3],
	float *b
) {
	*b=sqrtf(
		a[0]*a[0]+
		a[1]*a[1]+
		a[2]*a[2]
	);
}

static inline void vml_vec3_lerp(
	float a[3],
	float b[3],
	float t,
	float c[3]
) {
	float it = 1-t;

	c[0] = vml_lerp(a[0],b[0],t);
	c[1] = vml_lerp(a[1],b[1],t);
	c[2] = vml_lerp(a[2],b[2],t);
}

// QUATERNION FUNCTIONS

void vml_quat_lerp(
	float a[4],
	float b[4],
	float t,
	float c[4]
) {
	float ax = a[0];
	float ay = a[1];
	float az = a[2];
	float aw = a[3];
	float bx = b[0];
	float by = b[1];
	float bz = b[2];
	float bw = b[3];

	float am = sqrtf(ax*ax+ay*ay+az*az+aw*aw);
	float bm = sqrtf(bx*bx+by*by+bz*bz+bw*bw);

	ax /= am;
	ay /= am;
	az /= am;
	aw /= am;
	bx /= bm;
	by /= bm;
	bz /= bm;
	bw /= bm;

	float dot = ax*bx+ay*by+az*bz+aw*bw;

	if (dot<0) {
		bx  = -bx;
		by  = -by;
		bz  = -bz;
		bw  = -bw;
		dot = -dot;
	}

	if (dot>0.9995) {
		float cx = vml_lerp(ax,bx,t);
		float cy = vml_lerp(ay,by,t);
		float cz = vml_lerp(az,bz,t);
		float cw = vml_lerp(aw,bw,t);
		float cm = sqrtf(cx*cx+cy*cy+cz*cz+cw*cw);

		c[0] = cx/cm;
		c[1] = cy/cm;
		c[2] = cz/cm;
		c[3] = cw/cm;

		return;
	}

	float t0  = acosf(dot);
	float t1  = t0*t;
	float st0 = sin(t0);
	float st1 = sin(t1);
	float s0  = cosf(t1)-dot*st1/st0;
	float s1  = st1/st0;

	c[0] = ax*s0+bx*s1;
	c[1] = ay*s0+by*s1;
	c[2] = az*s0+bz*s1;
	c[3] = aw*s0+bw*s1;
}

// MAT3 FUNCTIONS

static inline void vml_mat3_id(
	float a[9]
) {
	a[0] = 1;
	a[1] = 0;
	a[2] = 0;
	a[3] = 0;
	a[4] = 1;
	a[5] = 0;
	a[6] = 0;
	a[7] = 0;
	a[8] = 1;
}

static inline void vml_mat3_mov(
	float a[9],
	float b[9]
) {
	for (unsigned int i=0; i<9; i++) {
		b[i] = a[i];
	}
}

void vml_mat3_transpose(
	float a[9],
	float b[9]
) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a10 = a[3];
	float a11 = a[4];
	float a12 = a[5];
	float a20 = a[6];
	float a21 = a[7];
	float a22 = a[8];

	b[0] = a00;
	b[1] = a10;
	b[2] = a20;
	b[3] = a01;
	b[4] = a11;
	b[5] = a21;
	b[6] = a02;
	b[7] = a12;
	b[8] = a22;
}

void vml_mat3_vec3_mul(
	float a[9],
	float b[3],
	float c[3]
) {
	float r0 = b[0]*a[0]+b[1]*a[1]+b[2]*a[2];
	float r1 = b[0]*a[3]+b[1]*a[4]+b[2]*a[5];
	float r2 = b[0]*a[6]+b[1]*a[7]+b[2]*a[8];

	c[0] = r0;
	c[1] = r1;
	c[2] = r2;
}

void vml_mat3_mul(
	float a[9],
	float b[9],
	float c[9]
) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a10 = a[3];
	float a11 = a[4];
	float a12 = a[5];
	float a20 = a[6];
	float a21 = a[7];
	float a22 = a[8];

	float b00 = b[0];
	float b01 = b[1];
	float b02 = b[2];
	float b10 = b[3];
	float b11 = b[4];
	float b12 = b[5];
	float b20 = b[6];
	float b21 = b[7];
	float b22 = b[8];

	c[0] = a00*b00+a01*b10+a02*b20;
	c[1] = a00*b01+a01*b11+a02*b21;
	c[2] = a00*b02+a01*b12+a02*b22;
	c[3] = a10*b00+a11*b10+a12*b20;
	c[4] = a10*b01+a11*b11+a12*b21;
	c[5] = a10*b02+a11*b12+a12*b22;
	c[6] = a20*b00+a21*b10+a22*b20;
	c[7] = a20*b01+a21*b11+a22*b21;
	c[8] = a20*b02+a21*b12+a22*b22;
}

void vml_mat3_inv(
	float a[9],
	float b[9]
) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a10 = a[3];
	float a11 = a[4];
	float a12 = a[5];
	float a20 = a[6];
	float a21 = a[7];
	float a22 = a[8];

	float det=(
		a00*(a11*a22-a21*a12)-
		a01*(a10*a22-a12*a20)+
		a02*(a10*a21-a11*a20)
	);

	b[0] = (a11*a22-a21*a12)/det;
	b[1] = (a02*a21-a01*a22)/det;
	b[2] = (a01*a12-a02*a11)/det;
	b[3] = (a12*a20-a10*a22)/det;
	b[4] = (a00*a22-a02*a20)/det;
	b[5] = (a10*a02-a00*a12)/det;
	b[6] = (a10*a21-a20*a11)/det;
	b[7] = (a20*a01-a00*a21)/det;
	b[8] = (a00*a11-a10*a01)/det;
}

void vml_mat3_scale_set(
	float a[9],
	float b[3],
	float c[9]
) {
	vml_mat3_mul(a,(float[9]){
		b[0], 0,    0,
		0,    b[1], 0,
		0,    0,    b[2]
	},c);
}

// MAT4 FUNCTIONS

static inline void vml_mat4_id(
	float a[16]
) {
	a[0]  = 1;
	a[1]  = 0;
	a[2]  = 0;
	a[3]  = 0;
	a[4]  = 0;
	a[5]  = 1;
	a[6]  = 0;
	a[7]  = 0;
	a[8]  = 0;
	a[9]  = 0;
	a[10] = 1;
	a[11] = 0;
	a[12] = 0;
	a[13] = 0;
	a[14] = 0;
	a[15] = 1;
}

static inline void vml_mat4_mov(
	float a[16],
	float b[16]
) {
	for (unsigned int i=0; i<16; i++) {
		b[i] = a[i];
	}
}


void vml_mat4_transpose(
	float a[16],
	float b[16]
) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a03 = a[3];
	float a10 = a[4];
	float a11 = a[5];
	float a12 = a[6];
	float a13 = a[7];
	float a20 = a[8];
	float a21 = a[9];
	float a22 = a[10];
	float a23 = a[11];
	float a30 = a[12];
	float a31 = a[13];
	float a32 = a[14];
	float a33 = a[15];

	b[0]  = a00;
	b[1]  = a10;
	b[2]  = a20;
	b[3]  = a30;
	b[4]  = a01;
	b[5]  = a11;
	b[6]  = a21;
	b[7]  = a31;
	b[8]  = a02;
	b[9]  = a12;
	b[10] = a22;
	b[11] = a32;
	b[12] = a03;
	b[13] = a13;
	b[14] = a23;
	b[15] = a33;
}

void vml_mat4_vec3_mul(
	float a[16],
	float b[3],
	float c[3]
) {
	float r0 = b[0]*a[0]+b[1]*a[1]+b[2]*a[2]+a[3];
	float r1 = b[0]*a[4]+b[1]*a[5]+b[2]*a[6]+a[7];
	float r2 = b[0]*a[8]+b[1]*a[9]+b[2]*a[10]+a[11];

	c[0] = r0;
	c[1] = r1;
	c[2] = r2;
}

void vml_mat4_mul(
	float a[16],
	float b[16],
	float c[16]
) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a03 = a[3];
	float a10 = a[4];
	float a11 = a[5];
	float a12 = a[6];
	float a13 = a[7];
	float a20 = a[8];
	float a21 = a[9];
	float a22 = a[10];
	float a23 = a[11];
	float a30 = a[12];
	float a31 = a[13];
	float a32 = a[14];
	float a33 = a[15];

	float b00 = b[0];
	float b01 = b[1];
	float b02 = b[2];
	float b03 = b[3];
	float b10 = b[4];
	float b11 = b[5];
	float b12 = b[6];
	float b13 = b[7];
	float b20 = b[8];
	float b21 = b[9];
	float b22 = b[10];
	float b23 = b[11];
	float b30 = b[12];
	float b31 = b[13];
	float b32 = b[14];
	float b33 = b[15];

	c[0]  = a00*b00+a01*b10+a02*b20+a03*b30;
	c[1]  = a00*b01+a01*b11+a02*b21+a03*b31;
	c[2]  = a00*b02+a01*b12+a02*b22+a03*b32;
	c[3]  = a00*b03+a01*b13+a02*b23+a03*b33;
	c[4]  = a10*b00+a11*b10+a12*b20+a13*b30;
	c[5]  = a10*b01+a11*b11+a12*b21+a13*b31;
	c[6]  = a10*b02+a11*b12+a12*b22+a13*b32;
	c[7]  = a10*b03+a11*b13+a12*b23+a13*b33;
	c[8]  = a20*b00+a21*b10+a22*b20+a23*b30;
	c[9]  = a20*b01+a21*b11+a22*b21+a23*b31;
	c[10] = a20*b02+a21*b12+a22*b22+a23*b32;
	c[11] = a20*b03+a21*b13+a22*b23+a23*b33;
	c[12] = a30*b00+a31*b10+a32*b20+a33*b30;
	c[13] = a30*b01+a31*b11+a32*b21+a33*b31;
	c[14] = a30*b02+a31*b12+a32*b22+a33*b32;
	c[15] = a30*b03+a31*b13+a32*b23+a33*b33;
}

void vml_mat4_inv(
	float a[16],
	float b[16]
) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a03 = a[3];
	float a10 = a[4];
	float a11 = a[5];
	float a12 = a[6];
	float a13 = a[7];
	float a20 = a[8];
	float a21 = a[9];
	float a22 = a[10];
	float a23 = a[11];
	float a30 = a[12];
	float a31 = a[13];
	float a32 = a[14];
	float a33 = a[15];

	float b00 =  a11*a22*a33-a11*a23*a32-a21*a12*a33+a21*a13*a32+a31*a12*a23-a31*a13*a22;
	float b01 = -a01*a22*a33+a01*a23*a32+a21*a02*a33-a21*a03*a32-a31*a02*a23+a31*a03*a22;
	float b02 =  a01*a12*a33-a01*a13*a32-a11*a02*a33+a11*a03*a32+a31*a02*a13-a31*a03*a12;
	float b03 = -a01*a12*a23+a01*a13*a22+a11*a02*a23-a11*a03*a22-a21*a02*a13+a21*a03*a12;
	float b10 = -a10*a22*a33+a10*a23*a32+a20*a12*a33-a20*a13*a32-a30*a12*a23+a30*a13*a22;
	float b11 =  a00*a22*a33-a00*a23*a32-a20*a02*a33+a20*a03*a32+a30*a02*a23-a30*a03*a22;
	float b12 = -a00*a12*a33+a00*a13*a32+a10*a02*a33-a10*a03*a32-a30*a02*a13+a30*a03*a12;
	float b13 =  a00*a12*a23-a00*a13*a22-a10*a02*a23+a10*a03*a22+a20*a02*a13-a20*a03*a12;
	float b20 =  a10*a21*a33-a10*a23*a31-a20*a11*a33+a20*a13*a31+a30*a11*a23-a30*a13*a21;
	float b21 = -a00*a21*a33+a00*a23*a31+a20*a01*a33-a20*a03*a31-a30*a01*a23+a30*a03*a21;
	float b22 =  a00*a11*a33-a00*a13*a31-a10*a01*a33+a10*a03*a31+a30*a01*a13-a30*a03*a11;
	float b23 = -a00*a11*a23+a00*a13*a21+a10*a01*a23-a10*a03*a21-a20*a01*a13+a20*a03*a11;
	float b30 = -a10*a21*a32+a10*a22*a31+a20*a11*a32-a20*a12*a31-a30*a11*a22+a30*a12*a21;
	float b31 =  a00*a21*a32-a00*a22*a31-a20*a01*a32+a20*a02*a31+a30*a01*a22-a30*a02*a21;
	float b32 = -a00*a11*a32+a00*a12*a31+a10*a01*a32-a10*a02*a31-a30*a01*a12+a30*a02*a11;
	float b33 =  a00*a11*a22-a00*a12*a21-a10*a01*a22+a10*a02*a21+a20*a01*a12-a20*a02*a11;

	float det = a00*b00+a01*b10+a02*b20+a03*b30;

	b[0]  = b00/det;
	b[1]  = b01/det;
	b[2]  = b02/det;
	b[3]  = b03/det;
	b[4]  = b10/det;
	b[5]  = b11/det;
	b[6]  = b12/det;
	b[7]  = b13/det;
	b[8]  = b20/det;
	b[9]  = b21/det;
	b[10] = b22/det;
	b[11] = b23/det;
	b[12] = b30/det;
	b[13] = b31/det;
	b[14] = b32/det;
	b[15] = b33/det;
}

void vml_mat4_mat3_get(
	float a[16],
	float b[9]
) {
	b[0] = a[0];
	b[1] = a[1];
	b[2] = a[2];
	b[3] = a[4];
	b[4] = a[5];
	b[5] = a[6];
	b[6] = a[8];
	b[7] = a[9];
	b[8] = a[10];
}

void vml_mat4_euler_set(
	float a[16],
	float b[3],
	float c[16]
) {
	float cx = cosf(a[0]);
	float cy = cosf(a[1]);
	float cz = cosf(a[2]);
	float sx = sinf(a[0]);
	float sy = sinf(a[1]);
	float sz = sinf(a[2]);

	b[0]  = cy*cz;
	b[1]  = -cy*sz;
	b[2]  = sy;
	b[3]  = 0;
	b[4]  = cz*sx*sy+cx*sz;
	b[5]  = cx*cz-sx*sy*sz;
	b[6]  = -cy*sx;
	b[7]  = 0;
	b[8]  = sx*sz-cx*cz*sy;
	b[9]  = cz*sx+cx*sy*sz;
	b[10] = cx*cy;
	b[11] = 0;
	b[12] = 0;
	b[13] = 0;
	b[14] = 0;
	b[15] = 1;
}

void vml_mat4_euler_get(
	float a[16],
	float b[3]
) {
	b[0] = atan2f(-a[6],a[10]);
	b[1] = asinf(a[2]);
	b[2] = atan2f(-a[1],a[0]);
}

void vml_mat4_quat_set(
	float a[16],
	float b[4],
	float c[16]
) {
	float x = b[0];
	float y = b[1];
	float z = b[2];
	float w = b[3];

	float xx = x*x;
	float yy = y*y;
	float zz = z*z;
	float xy = x*y;
	float xz = x*z;
	float yz = y*z;
	float xw = x*w;
	float yw = y*w;
	float zw = z*w;

	c[0]  = 1-2*yy-2*zz;
	c[1]  = 2*(xy-zw);
	c[2]  = 2*(xz+yw);
	c[3]  = a[3];
	c[4]  = 2*(xy+zw);
	c[5]  = 1-2*xx-2*zz;
	c[6]  = 2*(yz-xw);
	c[7]  = a[7];
	c[8]  = 2*(xz-yw);
	c[9]  = 2*(yz+xw);
	c[10] = 1-2*xx-2*yy;
	c[11] = a[11];
	c[12] = a[12];
	c[13] = a[13];
	c[14] = a[14];
	c[15] = a[15];
}

void vml_mat4_quat_get(
	float a[16],
	float b[4]
) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a10 = a[4];
	float a11 = a[5];
	float a12 = a[6];
	float a20 = a[8];
	float a21 = a[9];
	float a22 = a[10];

	float x = sqrtf(vml_max(( a00-a11-a22+1)/4,0));
	float y = sqrtf(vml_max((-a00+a11-a22+1)/4,0));
	float z = sqrtf(vml_max((-a00-a11+a22+1)/4,0));
	float w = sqrtf(vml_max(( a00+a11+a22+1)/4,0));

	if (w>=x && w>=y && w>=z) {
		x = x*vml_sign(a21-a12);
		y = y*vml_sign(a02-a20);
		z = z*vml_sign(a10-a01);
	} else if (x>=w && x>=y && x>=z) {
		y = y*vml_sign(a10+a01);
		z = z*vml_sign(a02+a20);
		w = w*vml_sign(a21-a12);
	} else if (y>=w && y>=x && y>=z) {
		x = x*vml_sign(a10+a01);
		z = z*vml_sign(a21+a12);
		w = w*vml_sign(a02-a20);
	} else if (z>=w && z>=x && z>=y) {
		x = x*vml_sign(a20+a02);
		y = y*vml_sign(a21+a12);
		w = w*vml_sign(a10-a01);
	}

	float m = sqrtf(x*x+y*y+z*z+w*w);

	b[0] = x/m;
	b[1] = y/m;
	b[2] = z/m;
	b[3] = w/m;
}

void vml_mat4_axis_set(
	float a[16],
	float b[4],
	float c[16]
) {
	float x = b[0];
	float y = b[1];
	float z = b[2];
	float w = b[3];

	float cw = cosf(w);
	float sw = sinf(w);
	float m  = sqrtf(x*x+y*y+z*z);

	x /= m;
	y /= m;
	z /= m;

	float xx = x*x;
	float yy = y*y;
	float zz = z*z;

	c[0]  = cw+xx*(1-cw);
	c[1]  = yy*(1-cw)-z*sw;
	c[2]  = zz*(1-cw)+y*sw;
	c[3]  = a[3];
	c[4]  = xx*(1-cw)+z*sw;
	c[5]  = cw+yy*(1-cw);
	c[6]  = zz*(1-cw)-x*sw;
	c[7]  = a[7];
	c[8]  = x*z*(1-cw)-y*sw;
	c[9]  = y*z*(1-cw)+x*sw;
	c[10] = cw+zz*(1-cw);
	c[11] = a[11];
	c[12] = a[12];
	c[13] = a[13];
	c[14] = a[14];
	c[15] = a[15];
}

void vml_mat4_axis_get(
	float a[16],
	float b[4]
) {
	float a21_a12 = a[9]-a[6];
	float a02_a20 = a[2]-a[8];
	float a10_a01 = a[4]-a[1];

	float m=sqrtf(
		a21_a12*a21_a12+
		a02_a20*a02_a20+
		a10_a01*a10_a01
	);

	b[0] = a21_a12/m;
	b[1] = a02_a20/m;
	b[2] = a10_a01/m;
	b[3] = acosf((a[0]+a[5]+a[10]-1)/2);
}

void vml_mat4_look_set(
	float a[16],
	float b[3],
	float c[3],
	float d[16]
) {
	float fx = b[0];
	float fy = b[1];
	float fz = b[2];
	float ux = c[0];
	float uy = c[1];
	float uz = c[2];
	float rx = uy*fz-uz*fy;
	float ry = uz*fx-ux*fz;
	float rz = ux*fy-uy*fx;
	float rm = sqrtf(rx*rx+ry*ry+rz*rz);

	rx /= rm;
	ry /= rm;
	rz /= rm;

	c[0]  = rx;
	c[1]  = fy*rz-fz*ry;
	c[2]  = fx;
	c[3]  = a[3];
	c[4]  = ry;
	c[5]  = fz*rx-fx*rz;
	c[6]  = fy;
	c[7]  = a[7];
	c[8]  = rz;
	c[9]  = fx*ry-fy*rx;
	c[10] = fz;
	c[11] = a[11];
	c[12] = a[12];
	c[13] = a[13];
	c[14] = a[14];
	c[15] = a[15];
}

void vml_mat4_look_get(
	float a[16],
	float b[3]
) {
	b[0] = a[2];
	b[1] = a[6];
	b[2] = a[10];
}

void vml_mat4_position_set(
	float a[16],
	float b[3],
	float c[16]
) {
	vml_mat4_mov(a,c);

	c[3]  = b[0];
	c[7]  = b[1];
	c[11] = b[2];
}

void vml_mat4_position_get(
	float a[16],
	float b[3]
) {
	b[0] = a[3];
	b[1] = a[7];
	b[2] = a[11];
}

void vml_mat4_scale_set(
	float a[16],
	float b[3],
	float c[16]
) {
	float rx = a[0];
	float ry = a[4];
	float rz = a[8];
	float tx = a[1];
	float ty = a[5];
	float tz = a[9];
	float fx = a[2];
	float fy = a[6];
	float fz = a[10];	

	float rm = b[0]*1/sqrtf(rx*rx+ry*ry+rz*rz);
	float tm = b[1]*1/sqrtf(tx*tx+ty*ty+tz*tz);
	float fm = b[2]*1/sqrtf(fx*fx+fy*fy+fz*fz);

	c[0]  = rx*rm;
	c[1]  = tx*tm;
	c[2]  = fx*fm;
	c[3]  = a[3];
	c[4]  = ry*rm;
	c[5]  = ty*tm;
	c[6]  = fy*fm;
	c[4]  = a[4];
	c[8]  = rz*rm;
	c[9]  = tz*tm;
	c[10] = fz*fm;
	c[11] = a[11];
	c[12] = a[12];
	c[13] = a[13];
	c[14] = a[14];
	c[15] = a[15];
}

void vml_mat4_scale_get(
	float a[16],
	float b[3]
) {
	float rx = a[0];
	float ry = a[4];
	float rz = a[8];
	float tx = a[1];
	float ty = a[5];
	float tz = a[9];
	float fx = a[2];
	float fy = a[6];
	float fz = a[10];	

	b[0] = sqrtf(rx*rx+ry*ry+rz*rz);
	b[1] = sqrtf(tx*tx+ty*ty+tz*tz);
	b[2] = sqrtf(fx*fx+fy*fy+fz*fz);
}

void vml_mat4_lerp(
	float a[16],
	float b[16],
	float t,
	float c[16]
) {
	float aq[4];
	float bq[4];
	float cq[4];
	float ap[3];
	float bp[3];
	float cp[3];

	vml_mat4_quat_get(a,aq);
	vml_mat4_quat_get(b,bq);

	vml_mat4_position_get(a,ap);
	vml_mat4_position_get(b,bp);

	vml_quat_lerp(aq,bq,t,cq);
	vml_vec3_lerp(ap,bp,t,cp);

	vml_mat4_quat_set(c,cq,c);
	vml_mat4_position_set(c,cp,c);
}

void vml_mat4_perspective(
	float v,
	float r,
	float n,
	float f,
	float a[16]
) {
	float t = tanf(v/2);

	a[0]  = 1/(t*r);
	a[1]  = 0;
	a[2]  = 0;
	a[3]  = 0;
	a[4]  = 0;
	a[5]  = 1/t;
	a[6]  = 0;
	a[7]  = 0;
	a[8]  = 0;
	a[9]  = 0;
	a[10] = -(f+n)/(f-n);
	a[11] = -(2*f*n)/(f-n);
	a[12] = 0;
	a[13] = 0;
	a[14] = -1;
	a[15] = 0;
}

void vml_mat4_orthographic(
	float l,
	float r,
	float t,
	float b,
	float n,
	float f,
	float a[16]
) {
	a[0] = 2/(r-l);
	a[1] = 0;
	a[2] = 0;
	a[3] = -(r+l)/(r-l);
	a[4] = 0;
	a[5] = 2/(t-b);
	a[6] = 0;
	a[7] = -(t+b)/(t-b);
	a[8] = 0;
	a[9] = 0;
	a[10] = -2/(f-n);
	a[11] = -(f+n)/(f-n);
	a[12] = 0;
	a[13] = 0;
	a[14] = 0;
	a[15] = 1;
}

#endif
#endif