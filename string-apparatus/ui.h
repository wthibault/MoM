//
// ui.h
//

#include <FTGL/ftgl.h>
#include <string>
class StringModel;

class UI {
 public:
  UI ( const char *fontpath );
  ~UI();

  virtual void printAt ( float x, float y, std::string msg );
  virtual void draw(int w, int h, StringModel *str);
  FTGLPixmapFont font;
};
