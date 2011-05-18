
OS_NAME := $(shell uname)

CXX = g++
CC = gcc
BASEDIR = .
INCLUDE = -I$(BASEDIR) $(shell wx-config --cxxflags)
OPTFLAGS =  -O3 -g3 -Wall -DNDEBUG
CXXFLAGS = $(OPTFLAGS) $(INCLUDE)

UTIL = Util/Util.cpp Util/ObjIO.cpp Util/ColorMap.cpp \
 Util/Heap.cpp Util/JetColorMap.cpp Util/HSVColorMap.cpp \
 Util/BlackWhiteColorMap.cpp Util/HotColorMap.cpp \
 Util/IsoContourColorMap.cpp Util/ColorMapFactory.cpp Util/trackball.cpp \
 Util/GLObjectPlayback.cpp

GUI =	main.cpp GUI/GUI.cpp GUI/FrameMain.cpp GUI/GLViewer.cpp \
			GUI/ScalarCutPlane.cpp GUI/VectorCutPlane.cpp \
			GUI/GLAxis.cpp  GUI/GLGridPlane.cpp GUI/GLObject.cpp \
			GUI/FluidVoxelCutPlane.cpp

MESH = Geometry/HalfEdgeMesh.cpp Decimation/DecimationMesh.cpp Decimation/SimpleDecimationMesh.cpp\
Decimation/QuadricDecimationMesh.cpp Util/MarchingCubes.cpp Geometry/SimpleMesh.cpp Geometry/Mesh.cpp \
Geometry/LineStrip.cpp

IMPLICITS =  Geometry/Implicit.cpp Geometry/Quadric.cpp Geometry/Sphere.cpp Geometry/SphereFractal.cpp \
 Geometry/SignedDistanceSphere.cpp Geometry/Cube.cpp Geometry/ImplicitMesh.cpp

LEVELSET = Levelset/LevelSetGrid.cpp Levelset/LevelSet.cpp Levelset/LevelSetOperator.cpp

FLUIDSOLVER = Fluid/FluidSolver.cpp

SUBDIVISION = Subdivision/LoopSubdivisionMesh.cpp \
							Subdivision/AdaptiveLoopSubdivisionMesh.cpp \
							Subdivision/UniformCubicSplineSubdivisionCurve.cpp \
							Subdivision/UniformCubicSpline.cpp


SOURCE =  $(UTIL) $(GUI) $(MESH) $(IMPLICITS) $(LEVELSET)\
$(FLUIDSOLVER) $(SUBDIVISION)

OBJ = $(SOURCE:.cpp=.o)

ifeq "$(OS_NAME)" "Darwin"
LDFLAGS = -framework GLUT -framework OpenGL -lobjc $(shell wx-config --libs --gl-libs)
else
LDFLAGS = -lGL -lGLU -lglut $(shell wx-config --libs --gl-libs)
endif

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o "main" $(OBJ)
	mkdir -p Moa.app
	mkdir -p Moa.app/Contents
	mkdir -p Moa.app/Contents/MacOS
	mkdir -p Moa.app/Contents/Resources
	cp Moa.plist Moa.app/Contents/Info.plist
	cp Moa.icns Moa.app/Contents/Resources
	echo -n 'APPL????' > Moa.app/Contents/PkgInfo
	cp main Moa.app/Contents/MacOS/Moa

# Automatic dependency updating
%.d: %.cpp
	@echo "Building dependencies for"  $(basename $@).o
	@$(CXX) -MM -MT $(basename $@).o $(INCLUDE) $< > $@

ifneq ($(MAKECMDGOALS),clean)
include $(SOURCE:.cpp=.d)
endif

clean:
	rm -fr *.o
	rm -fr *~
	rm -fr *.d
	rm -fr Decimation/*.o
	rm -fr Decimation/*~
	rm -fr Decimation/*.d
	rm -fr Fluid/*.o
	rm -fr Fluid/*~
	rm -fr Fluid/*.d
	rm -fr Geometry/*.o
	rm -fr Geometry/*~
	rm -fr Geometry/*.d
	rm -fr GUI/*.o
	rm -fr GUI/*~
	rm -fr GUI/*.d
	rm -fr Levelset/*.o
	rm -fr Levelset/*~
	rm -fr Levelset/*.d
	rm -fr Math/*.o
	rm -fr Math/*~
	rm -fr Math/*.d
	rm -fr Subdivision/*.o
	rm -fr Subdivision/*~
	rm -fr Subdivision/*.d
	rm -fr Util/*.o
	rm -fr Util/*~
	rm -fr Util/*.d
