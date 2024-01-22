#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CS248 {


// Implements SoftwareRenderer //

// pack color to uint 
unsigned int colorToUInt(const Color& c) {
        unsigned int r = 255 * c.r;
        unsigned int g = 255 * c.g;
        g <<= 8;
        unsigned int b = 255 * c.b;
        b <<= 16;
        unsigned int a = 255 * c.a;
        a <<= 24;
        return r | g | b | a;
}

// unpack uint from color
Color unitToColor(unsigned int uc) {
    float inv255 = 1.0 / 255.0;
    Color pixel_color;
    pixel_color.r = (uc & 0xFF) * inv255;
    pixel_color.g = ( (uc>>8) & 0xFF) * inv255;
    pixel_color.b = ((uc >> 16) & 0xFF) * inv255;
    pixel_color.a = ((uc >> 24) & 0xFF) * inv255;
    return pixel_color;
}

// fill a sample location with color
 void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color) {
  // Task 2: implement this function
  
  // check bounds
  if (sx < 0 || sx >= width*sample_rate) return;
  if (sy < 0 || sy >= height*sample_rate) return;

  Color pixel_color = unitToColor(sample_buffer[sx + sy * sbwidth()]);
  pixel_color = ref->alpha_blending_helper(pixel_color, color);
  sample_buffer[sx + sy * sbwidth()] = colorToUInt(pixel_color);
}

void SoftwareRendererImp::set_sample(int sx, int sy, const Color& color)
{
    unsigned int uc = colorToUInt(color);
    sample_buffer[sx + sy * sbwidth()] = uc;
}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// Task 2: Re-implement this function

	// check bounds
	if (x < 0 || x >= width) return;
	if (y < 0 || y >= height) return;

	
	pixel_buffer[4 * (x + y * width)] = (uint8_t) (color.r * 255);
	pixel_buffer[4 * (x + y * width) + 1] = (uint8_t)(color.g * 255);
	pixel_buffer[4 * (x + y * width) + 2] = (uint8_t)(color.b * 255);
	pixel_buffer[4 * (x + y * width) + 3] = (uint8_t)(color.a * 255);

}

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  // cout << transformation << "trans" << endl;
  //cout << "Starting draw svg" << endl;
  transformation = canvas_to_screen; // this is 3x3 matrix, initially transformation is identity 

  // canvas outline
  Vector2D a = transform(Vector2D(0, 0)); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0)); b.x++; b.y--;
  Vector2D c = transform(Vector2D(0, svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height)); d.x++; d.y++;

  svg_bbox_top_left = Vector2D(a.x+1, a.y+1);
  svg_bbox_bottom_right = Vector2D(d.x-1, d.y-1);

  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    transformation = Matrix3x3::identity(); // set transformation back to identity before drawing each element 
    draw_element(svg.elements[i]);
  }
  
  transformation = canvas_to_screen; 
  // draw canvas outline
  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to pixel buffer
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 2:
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  sample_buffer.resize(width * height * sample_rate * sample_rate);
  fill_sample_buffer(bgColor);
}

void SoftwareRendererImp::fill_sample_buffer(const Color& c)
{
    for (int x = 0; x < width * sample_rate; ++x) {
        for (int y = 0; y < height * sample_rate; ++y) {
            set_sample(x, y, c);
        }
    }
}

void SoftwareRendererImp::set_pixel_buffer( unsigned char* pixel_buffer,
                                             size_t width, size_t height ) {

  // Task 2:
  // You may want to modify this for supersampling support
  this->pixel_buffer = pixel_buffer;
  this->width = width;
  this->height = height;
  sample_buffer.resize(width * height * sample_rate * sample_rate);
  fill_sample_buffer(bgColor);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

	// Task 3 (part 1):
	// Modify this to implement the transformation stack
  // cout << element->transform << "element_transform" << endl;

  Matrix3x3 transformation0 = transformation; // original transfomration 
  Matrix3x3 transformation_g = transformation0.operator*(element->transform); // multiply transform together 
  transformation = canvas_to_screen.operator*(transformation_g); // respect to canvas_to_screen
  // cout << transformation << "trans" << endl;
  //cout <<"Drawing: "<< element->type << std::endl;
	switch (element->type) {
	case POINT:
		draw_point(static_cast<Point&>(*element));
		break;
	case LINE:
		draw_line(static_cast<Line&>(*element));
		break;
	case POLYLINE:
		draw_polyline(static_cast<Polyline&>(*element));
		break;
	case RECT:
		draw_rect(static_cast<Rect&>(*element));
		break;
	case POLYGON:
		draw_polygon(static_cast<Polygon&>(*element));
		break;
	case ELLIPSE:
		draw_ellipse(static_cast<Ellipse&>(*element));
		break;
	case IMAGE:
		draw_image(static_cast<Image&>(*element));
		break;
	case GROUP:
    transformation = transformation_g;
		draw_group(static_cast<Group&>(*element));  
		break;
	default:
		break;
	}
  transformation = transformation0; // set transform back to original transfomration 
}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) {

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;

  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));

  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Advanced Task
  // Implement ellipse rasterization

}

void SoftwareRendererImp::draw_image( Image& image ) {

  // Advanced Task
  // Render image element with rotation

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  // fill sample - NOT doing alpha blending!
  // TODO: Call fill_pixel here to run alpha blending
  pixel_buffer[4 * (sx + sy * width)] = (uint8_t)(color.r * 255);
  pixel_buffer[4 * (sx + sy * width) + 1] = (uint8_t)(color.g * 255);
  pixel_buffer[4 * (sx + sy * width) + 2] = (uint8_t)(color.b * 255);
  pixel_buffer[4 * (sx + sy * width) + 3] = (uint8_t)(color.a * 255);

}

void SoftwareRendererImp::bline(unsigned x1, unsigned y1,
          unsigned x2, unsigned y2, Color color)
{
  int dx  = x2 - x1,
          dy  = y2 - y1,
          y   = y1,
          eps = 0;

  float s = (float)dy / (float)dx; // calculate slope
  if ((dx == 0 && y1 < y2) || s > 0.0f) // positive slope
  {
    if (dx == 0 || s > 1.0f)
    {
      int x = x1;
      for ( int y = y1; y <= y2; y++ )  {
      fill_sample(x, y, color);
      eps += dx;
      if ( (eps << 1) >= dy )  {
        x++;  eps -= dy;
      }
      }
      return;
    } else 
    {
       for ( int x = x1; x <= x2; x++ )  {
        fill_sample(x, y, color);
        eps += dy;
        if ( (eps << 1) >= dx )  {
          y++;  eps -= dx;
        }
      }
      return;
    }
  }
  if ( s <0.0f) // negative slope
  {
    if (s>-1.0f)
    {
      float e = 0.0f;
      float y = y1;
      float m = s;
      for ( int x = x1; x <= x2; x++ )  {
        fill_sample(x, y, color);
        if (e+m > -0.5f)
        {
          e = e+m;
        }
        else
        {
          y-=1.0f;
          e = e + m +1.0f;
        }
      }
    } else {
      float e = 0.0f;
      float x = x1;
      float m = 1/s;
      for ( int y = y1; y >= y2; y-- )  {
        fill_sample(x, y, color);
        if (e+m > -0.5f)
        {
          e = e+m;
        }
        else
        {
          x+=1.0f;
          e = e + m +1.0f;
        }
      }
    }
    return;
  }

  for ( int x = x1; x <= x2; x++ )  {
    fill_sample(x, y, color);
    eps += dy;
    if ( (eps << 1) >= dx )  {
      y++;  eps -= dx;
    }
  }
}

Vector2D toSampleSpace(const Vector2D& pt, float sample_rate, const Vector2D& sc)
{
    Vector2D vc = (pt - sc) * sample_rate;

    return Vector2D(vc.x + sc.x * sample_rate, vc.y + sc.y * sample_rate);
}

Vector3D toSampleSpace(const Vector3D& vec, float sample_rate, const Vector2D& sc)
{
    Vector2D ss = toSampleSpace(Vector2D(vec.x, vec.y), sample_rate, sc);
    return Vector3D(ss.x, ss.y, 1.0f);
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 0:
  // Implement Bresenham's algorithm (delete the line below and implement your own)
  //ref->rasterize_line_helper(x0, y0, x1, y1, width, height, color, this);
  float dx = x1 - x0;
  float dy = y1 - y0;
  //std::cout<<"x0:"<<x0<<" y0:"<<y0<<" x1:"<<x1<<" y1:"<<y1<<"s: "<<dy/dx <<std::endl;

  if (x1 < x0) {
    swap(x1, x0);
    swap(y1, y0);
  }

  Vector2D sc((float)width / 2.0f, (float)height / 2.0f);
  Vector2D v0 = toSampleSpace(Vector2D(x0, y0), sample_rate, sc);
  Vector2D v1 = toSampleSpace(Vector2D(x1, y1), sample_rate, sc);
  bline(v0.x, v0.y, v1.x, v1.y, color);

  // Advanced Task
  // Drawing Smooth Lines with Line Width
}

void calcBoudingRect(const std::vector<Vector3D>& v, 
    Vector2D& min, Vector2D& max)
{
    min.x = v[0].x;
    min.y = v[0].y;
    max.x = v[0].x;
    max.y = v[0].y;
    for (int i = 1; i < v.size(); ++i) {
        if (v[i].x < min.x) {
            min.x = v[i].x;
        }
        if (v[i].x > max.x) {
            max.x = v[i].x;
        }
        if (v[i].y < min.y) {
            min.y = v[i].y;
        }
        if (v[i].y > max.y) {
            max.y = v[i].y;
        }
    }
    
}

void calcABC(const Vector3D& pt0,
    const Vector3D& pt1,
    const Vector3D& pt2, Vector3D& abc)
{
    abc.x = pt1.y - pt0.y;
    abc.y = -(pt1.x - pt0.x);
    abc.z = pt0.y * (pt1.x - pt0.x) - pt0.x * (pt1.y - pt0.y);
}

bool isInsideTri(float x, float y, const Vector3D& l0, const Vector3D& l1, const Vector3D& l2)
{
    Vector3D pxC(x, y, 1.0f);
    return dot(l0, pxC) <= 0.0f && dot(l1, pxC) <= 0.0f && dot(l2, pxC) <= 0.0f;
}

Vector2D toScreenSpace(const Vector2D& vec, float sample_rate)
{
    return Vector2D(vec.x / sample_rate, vec.y / sample_rate);
}

// rasterize a triangle implementation
void SoftwareRendererImp::drawTriImpl(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color)
{
    Vector2D sc((float)width / 2.0f, (float)height / 2.0f);
    Vector3D pt0 = toSampleSpace(Vector3D(x0 , y0 , 1.0f), sample_rate, sc);
    Vector3D pt1 = toSampleSpace(Vector3D(x1 , y1 , 1.0f), sample_rate, sc);
    Vector3D pt2 = toSampleSpace(Vector3D(x2 , y2 , 1.0f), sample_rate, sc);

    // Check is triangle is clockwise and flip it if necessary
    if (cross(Vector2D(x1 - x0, y1 - y0), Vector2D(x2 - x1, y2 - y1)) < 0)
    {
        swap(pt0, pt1);
    }
    
    //Calculate bounding rect
    Vector2D min, max;
    std::vector<Vector3D> v(3);
    v[0] = pt0;
    v[1] = pt1;
    v[2] = pt2;
    calcBoudingRect(v , min, max);
    Vector3D l0, l1, l2;
    
    // Calculate line coefficients
    calcABC(pt0, pt1, pt2, l0);
    calcABC(pt1, pt2, pt0, l1);
    calcABC(pt2, pt0, pt1, l2);
    //cout << "minX:" << (int)min.x << " minY:" << (int)min.y << " maxX:" << (int)max.x << " maxY:" << (int)max.y << endl;
    // Drawing triangle samples in the sample space 
    for (float x = (int)min.x; x < max.x; x += 1.0f)
    {
        for (float y = (int)min.y; y < max.y; y += 1.0f)
        {
            if (isInsideTri(x+0.5f, y+0.5f, l0, l1, l2))
            {
                fill_sample(x, y, color);
            }
        }
    }
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 1:
  // Implement triangle rasterization
  drawTriImpl(x0, y0, x1, y1, x2, y2, color);

  // Advanced Task
  // Implementing Triangle Edge Rules

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 4:
  // Implement image rasterization

}

// resolve samples to pixel buffer
void SoftwareRendererImp::resolve( void ) {

  // Task 2:
  // Implement supersampling
    for (int x = 0 /*svg_bbox_top_left.x*/; x < width/*svg_bbox_bottom_right.x*/; ++x)
    {
        for (int y = 0/*svg_bbox_top_left.y*/; y < height /*svg_bbox_bottom_right.y*/; ++y)
        {
            Color pc;
            for (int xS = x * sample_rate; xS < (x + 1) * sample_rate; ++xS)
            {
                for (int yS = y * sample_rate; yS < (y + 1) * sample_rate; ++yS)
                {
                    Color c = unitToColor(sample_buffer[xS + yS * sbwidth()]);
                    c *= 1.0f / (float)(sample_rate * sample_rate);
                    pc += c;
                }
            }
            fill_pixel(x, y, pc);
        }
    }

  fill_sample_buffer(bgColor);

  return;

}

Color SoftwareRendererImp::alpha_blending(Color pixel_color, Color color)
{
  // Task 5
  // Implement alpha compositing
  return pixel_color;
}


} // namespace CS248
