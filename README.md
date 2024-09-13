# Vector Matrix Library
![LICENSE](https://img.shields.io/badge/LICENSE-MIT-green.svg)

VML is a simple 3D math library intended for game development.

## Features
- Single header library written in C99
- 2D and 3D vector functions
- 3x3 and 4x4 matrix functions

## How to use
Include library
```c
#define VML_IMPLEMENTATION
#include "vml.h"
```
Adding two vec3's
```c
float vector_a[3] = {1,2,3};
float vector_b[3] = {4,5,6};

vml_vec3_add(vector_a,vector_b,vector_a); // vector_a = vector_a+vector_b
```
Multiplying mat4 with vec3
```c
float matrix_a[16] = VML_MAT4_IDENTITY;
float vector_b[3]  = {1,2,3};

vml_mat4_vec3_mul(matrix_a,vector_b,vector_b}; // vector_b = matrix_a*vector_b
```
Refer to [vml.h](vml.h) for to see rest of the functions.

## License
This software is free to use. You can modify it and redistribute it under the terms of the 
MIT license. Check [LICENSE](LICENSE) for further details.
