#ifndef STRINGMODELPRIMITIVE_H
#define STRINGMODELPRIMITIVE_H
//////
/////////////////////////////////////////////////////////////////////////
// StringModelPrimitive - an SSG rendering Primitive subclassed from ParticleSystem
//

class StringModelPrimitive : public ParticleSystem
{
public:
  StringModelPrimitive( StringModel *stringModel ) 
    : theString_ ( stringModel ), renderScale_ ( 10.0 )
  {
    StringModelPrimitive::init();
  }
  ~StringModelPrimitive () {}

  void init() 
  {
    points_.clear();
    normals_.clear();
    indices_.clear();
    texCoords_.clear();
    for ( int i = 0; i < theString_->numMasses; i++ ) {
      points_.push_back ( glm::vec3(0,0,0) );
      normals_.push_back ( glm::vec3(0,0,1) );
      indices_.push_back ( i );
      texCoords_.push_back ( glm::vec2(i/theString_->numMasses, 0) );
    }
    drawingPrimitive_ = GL_POINTS;
    Primitive::init();
  }

  void update(float dt) 
  {
  }

  void draw ( glm::mat4 mv, glm::mat4 proj, Material *mat ) 
  {
    // copy the positions of masses to the vertex array
    glBindVertexArray ( vao_ );

    points_.clear();
    indices_.clear();
    normals_.clear();
    texCoords_.clear();

    float deltaX = 1.0 / theString_->numMasses;
    for ( int i = 0; i < theString_->numMasses; i++ ) {
      float x,y,z;
      x = 2*(i * deltaX)-1;
      y = renderScale_ * theString_->yold[i];
      z = 0;
      points_.push_back( glm::vec3 ( x,y,z ) );
      normals_.push_back ( glm::vec3(0,0,1) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 0) );
      indices_.push_back ( i );
    }
    
    long int sizeofPoints = sizeof(glm::vec3)*points_.size();
    int sizeofNormals = sizeof(glm::vec3)*normals_.size();
    int sizeofTexCoords = sizeof(glm::vec2)*texCoords_.size();
    drawingPrimitive_ = GL_POINTS;

    glBindBuffer ( GL_ARRAY_BUFFER, arrayBuffer_ );
    glBufferSubData( GL_ARRAY_BUFFER, 0, sizeofPoints, &points_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints, sizeofNormals, &normals_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints + sizeofNormals, sizeofTexCoords, &texCoords_[0] );
    glBindBuffer ( GL_ARRAY_BUFFER, 0 );
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBuffer_);
    int sizeofIndices = indices_.size()*sizeof(unsigned int);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeofIndices, &indices_[0]);
    glBindBuffer ( GL_ELEMENT_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);
    
    Primitive::draw ( mv, proj, mat );

    //    printParams();
  }

  void printParam ( const std::string name, float value ) {
    std::cout << name << std::setw( 6 ) << value << std::endl;
  }

  void printParams() {
    // temp until HUD
    std::cout << std::endl << std::endl << std::endl;
    printParam ("Ktension    ", theString_->Ktension );
#ifdef NEW_STRING_MODEL
    printParam ("massDensity ", theString_->massDensity );
    printParam ("decayTime   ", theString_->decayTime );
#else
    printParam ("Kdamping    ", theString_->Kdamping );
#endif
    printParam ("vib freq.   ", theString_->vibratorFreq );
    printParam ("vib amp.    ", theString_->vibratorAmplitude );
    
    std::cout << std::endl;
    std::cout << "t/g tension (coarse) " << std::endl;
    std::cout << "T/G tension (fine) " << std::endl;
    std::cout << "e/d damping (coarse) " << std::endl;
    std::cout << "E/D damping (fine) " << std::endl;
    std::cout << "u/j vib freq (coarse) " << std::endl;
    std::cout << "U/J vib freq (fine) " << std::endl;
    std::cout << "i/k vib amp (coarse) " << std::endl;
    std::cout << "I/K vib amp (fine) " << std::endl;
    std::cout << "p/P pluck" << std::endl;
    std::cout << "v toggle vibrator " << std::endl;
    std::cout << "r reset" << std::endl;
    std::cout << "ESC exit" << std::endl;
  }

  StringModel *theString_;
  float        renderScale_;
};


/////////////////////////////////////////////////////////////////

class StringModelHistogramPrimitive : public ParticleSystem
{
public:
  StringModelHistogramPrimitive( StringModel *stringModel ) 
    : theString_ ( stringModel ), renderScale_ ( 8.0 )
  {
    StringModelHistogramPrimitive::init();
  }
  ~StringModelHistogramPrimitive () {}

  void init() 
  {
    points_.clear();
    normals_.clear();
    indices_.clear();
    texCoords_.clear();
    for ( int i = 0; i < theString_->numMasses; i++ ) {
      points_.push_back ( glm::vec3(0,0,0) );
      points_.push_back ( glm::vec3(0,0,0) );
      normals_.push_back ( glm::vec3(0,0,1) );
      normals_.push_back ( glm::vec3(0,0,1) );
      indices_.push_back ( i );
      indices_.push_back ( i );
      texCoords_.push_back ( glm::vec2(i/theString_->numMasses, 0) );
      texCoords_.push_back ( glm::vec2(i/theString_->numMasses, 1) );
    }
    drawingPrimitive_ = GL_LINES;
    Primitive::init();
  }

  void update(float dt) 
  {
  }

  void draw ( glm::mat4 mv, glm::mat4 proj, Material *mat ) 
  {
    // draw a line based on histogram of positions
    glBindVertexArray ( vao_ );

    points_.clear();
    indices_.clear();
    normals_.clear();
    texCoords_.clear();

    float deltaX = 1.0 / theString_->numMasses;
    for ( int i = 0; i < theString_->numMasses; i++ ) {
      float x,y0,y1,z;
      x = 2*(i * deltaX)-1;
      y0 = renderScale_ * theString_->histograms[i].minVal;
      y1 = renderScale_ * theString_->histograms[i].maxVal;
      z = 0;
      points_.push_back( glm::vec3 ( x,y0,z ) );
      points_.push_back( glm::vec3 ( x,y1,z ) );
      normals_.push_back ( glm::vec3(0,0,1) );
      normals_.push_back ( glm::vec3(0,0,1) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 0) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 1) );
      indices_.push_back ( 2*i );
      indices_.push_back ( 2*i+1 );
      theString_->histograms[i].clear();
      //      theString_->histograms[i].minVal = -1e-6;
      //      theString_->histograms[i].maxVal = +1e-6;
      theString_->histograms[i].minVal = 100;
      theString_->histograms[i].maxVal = -100;
    }

    
    long int sizeofPoints = sizeof(glm::vec3)*points_.size();
    int sizeofNormals = sizeof(glm::vec3)*normals_.size();
    int sizeofTexCoords = sizeof(glm::vec2)*texCoords_.size();
    drawingPrimitive_ = GL_LINES;

    glBindBuffer ( GL_ARRAY_BUFFER, arrayBuffer_ );
    glBufferSubData( GL_ARRAY_BUFFER, 0, sizeofPoints, &points_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints, sizeofNormals, &normals_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints + sizeofNormals, sizeofTexCoords, &texCoords_[0] );
    glBindBuffer ( GL_ARRAY_BUFFER, 0 );
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBuffer_);
    int sizeofIndices = indices_.size()*sizeof(unsigned int);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeofIndices, &indices_[0]);
    glBindBuffer ( GL_ELEMENT_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);
    
    Primitive::draw ( mv, proj, mat );

    //    printParams();
  }

  void printParam ( const std::string name, float value ) {
    std::cout << name << std::setw( 6 ) << value << std::endl;
  }

  void printParams() {
    // temp until HUD
    std::cout << std::endl << std::endl << std::endl;
    printParam ("Ktension    ", theString_->Ktension );
#ifndef NEW_STRING_MODEL
    printParam ("Kdamping    ", theString_->Kdamping );
#endif
    printParam ("vib freq.   ", theString_->vibratorFreq );
    printParam ("vib amp.    ", theString_->vibratorAmplitude );
    
    std::cout << std::endl;
    std::cout << "t/g tension (coarse) " << std::endl;
    std::cout << "T/G tension (fine) " << std::endl;
    std::cout << "e/d damping (coarse) " << std::endl;
    std::cout << "E/D damping (fine) " << std::endl;
    std::cout << "u/j vib freq (coarse) " << std::endl;
    std::cout << "U/J vib freq (fine) " << std::endl;
    std::cout << "i/k vib amp (coarse) " << std::endl;
    std::cout << "I/K vib amp (fine) " << std::endl;
    std::cout << "p/P pluck" << std::endl;
    std::cout << "v toggle vibrator " << std::endl;
    std::cout << "r reset" << std::endl;
    std::cout << "ESC exit" << std::endl;
  }

  StringModel *theString_;
  float        renderScale_;
};

#endif
