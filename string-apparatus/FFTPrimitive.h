#ifndef FFTPRIMITIVE_H
#define FFTPRIMITIVE_H

class FFTPrimitive : public ParticleSystem
{
public:
  FFTPrimitive( StringModel *stringModel ) 
    : theString_ ( stringModel )
  {
    FFTPrimitive::init();
  }
  ~FFTPrimitive () {}

  void init() 
  {
    points_.clear();
    normals_.clear();
    indices_.clear();
    texCoords_.clear();
    int nFreqBins = theString_->numFramesToAnalyze / 2;
    float step = 1.0 / nFreqBins;
    for ( int i = 0; i < nFreqBins; i++ ) {
      points_.push_back ( glm::vec3(i*step, 0.0, 0) );
      points_.push_back ( glm::vec3(i*step, 1.0, 0) );
      normals_.push_back ( glm::vec3(0,0,1) );
      normals_.push_back ( glm::vec3(0,0,1) );
      indices_.push_back ( 2*i );
      indices_.push_back ( 2*i+1 );
      texCoords_.push_back ( glm::vec2(i*step, 0) );
      texCoords_.push_back ( glm::vec2(i*step, 1) );
    }
    drawingPrimitive_ = GL_LINES;
    Primitive::init();
  }

  void update(float dt) 
  {
    // setup the fft
    // XXX move out of string into here XXX
    // copy them backwards, same diff to the fft.
    theString_->ringBuffer.copyReversed ( theString_->fftwIn, theString_->numFramesToAnalyze, 2 );

    // run the fft on them
    fftw_execute ( theString_->fftwPlan );

    // load the geometry into the vao
    points_.clear();
    indices_.clear();
    normals_.clear();
    texCoords_.clear();

    int nFreqBins = theString_->numFramesToAnalyze / 2;
    float deltaX = 1.0 / nFreqBins;
    for ( int i = 0; i < nFreqBins; i++ ) {
      float x,y,z;
      x = i * deltaX;

      //#define FFT_AMP_IN_DECIBELS
#ifdef FFT_AMP_IN_DECIBELS
      //      y = 20*log10(max(1e-4,fabs(theString_->fftwOut[i][0])));
      y = 20*log10(theString_->fftwOut[i][0]);
      y += 50;
#else
      y = abs(theString_->fftwOut[i][0]);
#endif

      z = 0;
      points_.push_back( glm::vec3 ( x,0,z ) );
      points_.push_back( glm::vec3 ( x,y,z ) );
      normals_.push_back ( glm::vec3(0,0,1) );
      normals_.push_back ( glm::vec3(0,0,1) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 0) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 1) );
      indices_.push_back ( 2*i );
      indices_.push_back ( 2*i+1 );
    }
  }

  void draw ( glm::mat4 mv, glm::mat4 proj, Material *mat ) 
  {
    // copy the positions of masses to the vertex array
    glBindVertexArray ( vao_ );

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
    
    //    glDisable ( GL_DEPTH_TEST );
    glLineWidth(2.0);
    Primitive::draw ( mv, proj, mat );
    //    glEnable ( GL_DEPTH_TEST );

    //    printParams();
  }

  StringModel *theString_;
  
};

#endif
