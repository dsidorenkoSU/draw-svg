#include "viewport.h"

#include "CS248.h"
#include <iostream>

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 3 (part 2): 
  // Set svg to normalized device coordinate transformation. 
  // Your input arguments are defined as SVG canvans coordinates.

  this->x = x;
  this->y = y;
  this->span = span; 

  double data[9] = {span*2, 0, x-span, 0, span*2, y-span, 0, 0, 1}; // construct matrix; this transform normalized space to x, y
  Matrix3x3 transform_matrix(data); // initialize matrix3x3 
  set_canvas_to_norm(transform_matrix.inv()); // inverse the matrix so it transforms x, y to normalized space 
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248
