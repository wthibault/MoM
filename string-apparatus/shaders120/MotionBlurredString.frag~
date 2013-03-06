#version 120
//
// MotionBlurredString.frag
//
uniform sampler2D Texture;
uniform vec4 AmbientProduct;
uniform vec4 DiffuseProduct;
uniform vec4 SpecularProduct;
uniform float Shininess;

// these are all in eye coords
varying vec3 Light, View, Normal;


const int  numBins = 256;
uniform float bins[numBins];
varying float u;


void main() 
{ 
  vec4 color;
  vec3 L = normalize(Light);
  vec3 V = normalize(View);
  vec3 N = normalize(Normal);
  vec3 H = normalize ( L + V );

  // Compute terms in the illumination equation
  vec4 ambient = AmbientProduct;
  
  // two-sided lighting
  float ldotn = dot ( L, N );
  float Kd = max( ldotn, -ldotn );
  vec4  diffuse = Kd * DiffuseProduct;
  
  float Ks = pow( max(dot(N, H), -dot(N,H)), Shininess );
  vec4  specular = Ks * SpecularProduct;
  
  int bin = (numBins-1) * u;
  gl_FragColor =  diffuse +  (ambient + specular);
  gl_FragColor.a = bins[bin];
} 
