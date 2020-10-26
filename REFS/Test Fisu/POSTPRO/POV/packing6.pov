#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"
#include "screen.inc"

global_settings {
  max_trace_level 64
}

camera {
  location < 0.03500000, 0.03150000,-0.04846653>
  look_at < 0.00000000, 0.01346653, 0.00000000>
  rotate<0,0,0>
}

background{rgb 1}

light_source{<-8,30,-20> color rgb <1,1,1>}
light_source{<25,12,-12> color rgb <1,1,1>}

#declare fwstd=finish{diffuse 0.4 specular 0.8 phong 0.8 reflection 0.03}
#declare pwstd=pigment{color rgbf<0.85,0.85,0.85,0.9>}
#declare fpstd=finish{reflection 0.1 specular 0.3 ambient 0.42}
#declare ppstd=pigment{rgbt <0.45,0.45,0.45, 0.60000000>}
#declare f1=finish{reflection 0.17 specular 0.3 ambient 0.42}
#declare t0=texture{pigment{rgbt <0.6,0.6,0.6, 0.60000000>} finish{f1}}
#declare t1=texture{pigment{rgbt <0.7,0.3,1.0, 0.60000000>} finish{f1}}
#declare t2=texture{pigment{rgbt <0.3,0.7,1.0, 0.60000000>} finish{f1}}
#declare t3=texture{pigment{rgbt <0.2,0.9,0.9, 0.60000000>} finish{f1}}
#declare t4=texture{pigment{rgbt <0.3,1.0,0.7, 0.60000000>} finish{f1}}
#declare t5=texture{pigment{rgbt <0.7,1.0,0.3, 0.60000000>} finish{f1}}
#declare t6=texture{pigment{rgbt <0.9,0.9,0.2, 0.60000000>} finish{f1}}
#declare t7=texture{pigment{rgbt <1.0,0.7,0.3, 0.60000000>} finish{f1}}
#declare t8=texture{pigment{rgbt <1.0,0.3,0.7, 0.60000000>} finish{f1}}
#declare t9=texture{pigment{rgbt <0.9,0.2,0.9, 0.60000000>} finish{f1}}
#declare t10=texture{pigment{rgbt <1.0,1.0,1.0, 0.60000000>} finish{f1}}


union{
#include "./POSTPRO/POV/particles6.pov"
  rotate <-90,0,0>
  pigment{ppstd} finish{fpstd}
}

union{
#include "./POSTPRO/POV/walls6.pov"
  rotate <-90,0,0>
  texture{
    pigment{pwstd}
    finish{fwstd}
  }
}
union{
#include "./POSTPRO/POV/borders6.pov"
  rotate <-90,0,0>
  texture{
    pigment { rgbt <0.000000,0.000000,0.000000,0.000000>}
    finish{ambient 0.6}
  }
}
