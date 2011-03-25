
#include "GLObject.h"
#include <iostream>

void GLObject::Render() {
	GLenum error = glGetError();
	if (error != GL_NO_ERROR) {
		const GLubyte * errorString = gluErrorString(error);
		std::cerr << "OpenGL Error: " << errorString << std::endl;
	}
}
