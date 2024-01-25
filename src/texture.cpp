#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Advanced Task
  // Implement mipmap for trilinear filtering

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 4: Implement nearest neighbour interpolation
  
  int tex_y = int(round(v*(tex.height)-0.5));  // get nearest in y
  int tex_x = int(round(u*(tex.width)-0.5));   // get nearest in x 

  // check boundary 
  if (tex_x > tex.width-1) {
    tex_x = tex.width-1; 
  }
  if (tex_y > tex.height-1) {
    tex_y = tex.height-1; 
  }
  if (tex_x < 0) { // < 0 case shouldnt happen...
    tex_x = 0; 
  }
  if (tex_y < 0) { // < 0 case shouldnt happen...
    tex_y = 0; 
  }
  
  // // check invalid level and return mangeta. this shouldnt happen...
  // if (tex_x > tex.width || tex_y > tex.width || tex_x < 0 || tex_y < 0) {
  //   return Color(1,0,1,1);
  // }
  int loc = (tex_y*tex.width+tex_x)*4; // row major so do y*width, *4 because of RGBA 

  // Divide 255.f because texture map store 0-255 but Color takes 0-1...no gamma correction 
  Color c(float(tex.mipmap[level].texels[loc])/255.f, float(tex.mipmap[level].texels[loc+1])/255.f, float(tex.mipmap[level].texels[loc+2])/255.f, float(tex.mipmap[level].texels[loc+3])/255.f); 
  // cout << c << endl;

  return c;
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 4: Implement bilinear filtering  
  float tex_y = v*float(tex.height)-0.5; // -1 because it runs from 0 to height-1 
  float tex_x = u*float(tex.width)-0.5;  

  // // check invalid level and return mangeta 
  // if (tex_x > tex.width || tex_y > tex.width || tex_x < 0 || tex_y < 0) {
  //   cout << u << " " << v << endl; 

  //   return Color(1,0,1,1);
  // }
  
  // check boundary 
  if (tex_x > tex.width-1) {
    tex_x = tex.width-1; 
  }
  if (tex_y > tex.height-1) {
    tex_y = tex.height-1; 
  }
  if (tex_x < 0) { // < 0 case shouldnt happen...
    tex_x = 0; 
  }
  if (tex_y < 0) { // < 0 case shouldnt happen...
    tex_y = 0; 
  }

  // calculate integer location between tex_x, text_y 
  float x0 = floor(tex_x); 
  float x1 = floor(tex_x)+1; 
  x1 = (x1 > tex.width) ? x0 : x1;  
  float y0 = floor(tex_y); 
  float y1 = floor(tex_y)+1; 
  y1 = (y1 > tex.height) ? y0 : y1;  

  // calculate location in the texture 
  int loc_x0y0 = (int(y0)*tex.width+int(x0))*4; // row major so do y*width, *4 because of RGBA 
  int loc_x1y0 = (int(y0)*tex.width+int(x1))*4; 
  int loc_x0y1 = (int(y1)*tex.width+int(x0))*4; 
  int loc_x1y1 = (int(y1)*tex.width+int(x1))*4; 

  // bilinear interpolation 
  float color_interp[4] = {0, 0, 0, 0}; 
  for (int i = 0; i < 4; i++ ) { //RGBA
    // get color  
    float color_x0y0 = float(tex.mipmap[level].texels[loc_x0y0+i]);
    float color_x1y0 = float(tex.mipmap[level].texels[loc_x1y0+i]);
    float color_x0y1 = float(tex.mipmap[level].texels[loc_x0y1+i]);
    float color_x1y1 = float(tex.mipmap[level].texels[loc_x1y1+i]);

    // interpolate in x first 
    float color_temp_y0 = (tex_x-x0)/(x1-x0)*color_x1y0 + (x1-tex_x)/(x1-x0)*color_x0y0; 
    float color_temp_y1 = (tex_x-x0)/(x1-x0)*color_x1y1 + (x1-tex_x)/(x1-x0)*color_x0y1; 
    // interpolate in y 
    float color_temp = (tex_y-y0)/(y1-y0)*(color_temp_y1) + (y1-tex_y)/(y1-y0)*(color_temp_y0); 

    // Divide 255.f because texture map store 0-255 but Color takes 0-1...no gamma correction 
    color_interp[i] = color_temp/255.f; 
  }  
  // cout << color_interp[0] << " " << color_interp[1] << " " << color_interp[2] << " " << color_interp[3] << endl;
  return Color(color_interp[0], color_interp[1], color_interp[2], color_interp[3]); 

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Advanced Task
  // Implement trilinear filtering

  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CS248
