#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"
#include "screen.inc"

global_settings {
  max_trace_level 64
}

camera {
  location <20.00000000,20.00000000,21.75000000>
  right 0.24*x*image_width/image_height
  up 0.25*y
  look_at < 0.00000000, 1.75000000, 0.00000000>
  rotate<90,0,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <1,1,1>}
light_source{<25,12,-12> color rgb <1,1,1>}

#declare f1=finish{reflection 0.17 specular 0.3 ambient 0.42}
#declare t0=texture{pigment{rgb <0.6,0.6,0.6>} finish{f1}}
#declare t1=texture{pigment{rgb <0.7,0.3,1>} finish{f1}}
#declare t2=texture{pigment{rgb <0.3,0.7,1>} finish{f1}}
#declare t3=texture{pigment{rgb <0.2,0.9,0.9>} finish{f1}}
#declare t4=texture{pigment{rgb <0.3,1,0.7>} finish{f1}}
#declare t5=texture{pigment{rgb <0.7,1,0.3>} finish{f1}}
#declare t6=texture{pigment{rgb <0.9,0.9,0.2>} finish{f1}}
#declare t7=texture{pigment{rgb <1,0.7,0.3>} finish{f1}}
#declare t8=texture{pigment{rgb <1,0.3,0.7>} finish{f1}}
#declare t9=texture{pigment{rgb <0.9,0.2,0.9>} finish{f1}}
#declare t10=texture{pigment{rgb <1,1,1>} finish{f1}}


union{
#include "particles100.pov"
  rotate <0,0,0>
  pigment{rgb 0.45} finish{reflection 0.1 specular 0.3 ambient 0.42}
}

union{
#include "walls100.pov"
  rotate <0,0,0>
  texture{
    pigment{color rgbf<0.8,0.8,0.8,0.9>}
    finish{diffuse 0.2 specular 0.8 phong 0.8 reflection 0.1}
  }
}
