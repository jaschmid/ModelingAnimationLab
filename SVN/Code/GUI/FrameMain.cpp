#include "FrameMain.h"

FrameMain::FrameMain(wxWindow* parent) : BaseFrameMain(parent)
{
  GLGridPlane * plane = new GLGridPlane("Grid");
  plane->SetDimensions(4.0, 4.0);
  AddUniqueObject(plane);

  //GLAxis * axis = new GLAxis("axis");
  //AddUniqueObject(axis);

  // Connect the "object selected" event triggered by GLViewer
  Connect(wxID_ANY, wxEVT_GL_OBJECT_SELECTED, wxCommandEventHandler(FrameMain::ObjectSelected));

  // Add the panels to show for all corresponding types of GLObjects
  //Lab1->
  mPanelSwitches[typeid(Mesh).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(Mesh).name()].push_back(mPanelVisualization);
  //Lab1<-
  //Lab2->
  mPanelSwitches[typeid(DecimationMesh).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(DecimationMesh).name()].push_back(mPanelVisualization);
  mPanelSwitches[typeid(DecimationMesh).name()].push_back(mPanelDecimation);
  //Lab2<-
  //Lab3->
  mPanelSwitches[typeid(UniformCubicSpline).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(UniformCubicSplineSubdivisionCurve).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(UniformCubicSplineSubdivisionCurve).name()].push_back(mPanelSubdivision);
  mPanelSwitches[typeid(LoopSubdivisionMesh).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(LoopSubdivisionMesh).name()].push_back(mPanelVisualization);
  mPanelSwitches[typeid(LoopSubdivisionMesh).name()].push_back(mPanelSubdivision);
  //Lab3<-
  //Lab4->
  mPanelSwitches[typeid(Implicit).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(Implicit).name()].push_back(mPanelVisualization);
  mPanelSwitches[typeid(Implicit).name()].push_back(mPanelImplicit);
  mPanelSwitches[typeid(ScalarCutPlane).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(ScalarCutPlane).name()].push_back(mPanelVisualization);
  mPanelSwitches[typeid(VectorCutPlane).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(VectorCutPlane).name()].push_back(mPanelVisualization);
  //Lab4<-

  // Add all available color maps
  std::list<std::string> maps = ColorMapFactory::GetColorMaps();
  std::list<std::string>::iterator iter = maps.begin();
  std::list<std::string>::iterator iend = maps.end();
  while (iter != iend) {
    mColorMapChoice->Append(wxString((*iter).c_str(), wxConvUTF8));
    iter++;
  }

  UpdatePanels();
  mGLViewer->Render();
}


void FrameMain::AddUniqueObject(GLObject * object)
{
  unsigned int i = 1;
  std::string name = object->GetName();
  while (!mGLViewer->AddObject(object)) {
    i++;
    std::stringstream s;
		s << i;
    object->SetName(name + " (" + s.str() + ")");
  }
  mObjectList->Append(wxString(object->GetName().c_str(), wxConvUTF8));
}


void FrameMain::RemoveObject(GLObject * object)
{
  mGLViewer->RemoveObject(object->GetName());
  mObjectList->Delete(mObjectList->FindString(wxString(object->GetName().c_str(), wxConvUTF8)));
}


void FrameMain::DeleteObjects( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  while (!objects.empty()) {
    GLObject * object = objects.back();
    RemoveObject(object);
    DeleteDependentObjects(object);
    delete object;
    objects.pop_back();
  }
  UpdatePanels();
  mGLViewer->Render();
}


void FrameMain::DeleteDependentObjects(GLObject * object)
{
  std::list<std::string> & objects = mDependentObjects[object->GetName()];
  std::list<std::string>::const_iterator iter = objects.begin();
  std::list<std::string>::const_iterator iend = objects.end();
  while (iter != iend) {
    GLObject * dependent = mGLViewer->GetObject(*iter);
    if (dependent != NULL) {
      RemoveObject(dependent);
      delete dependent;
    }
    iter++;
  }
}


void FrameMain::UpdateDependentObjects(GLObject * object)
{
  std::list<std::string> & objects = mDependentObjects[object->GetName()];
  std::list<std::string>::const_iterator iter = objects.begin();
  std::list<std::string>::const_iterator iend = objects.end();
  while (iter != iend) {
    Geometry * dependent = dynamic_cast<Geometry *>(mGLViewer->GetObject(*iter));
    if (dependent != NULL)
      dependent->Update();
    iter++;
  }
}


void FrameMain::ObjectSelected(wxCommandEvent & event)
{
  //mObjectList->SetSelection(wxNOT_FOUND);
	for(size_t i = 0; i < mObjectList->GetCount(); i++)
    mObjectList->Deselect(i);

  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    mObjectList->SetStringSelection(wxString((*iter)->GetName().c_str(), wxConvUTF8));
    iter++;
  }

  UpdatePanels();
}


void FrameMain::SelectObjects( wxCommandEvent& event )
{
  mGLViewer->DeselectAllObjects();
  wxArrayInt objects;
  int numObjects = mObjectList->GetSelections(objects);
  for (int i = 0; i < numObjects; i++)
    mGLViewer->SelectObject(std::string(mObjectList->GetString(objects[i]).mb_str()));
}


void FrameMain::MoveObjectsUp( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    mGLViewer->MoveObject((*iter)->GetName(), -1);

    wxString name((*iter)->GetName().c_str(), wxConvUTF8);
    int i = mObjectList->FindString(name);
    if (i > 0) {
      mObjectList->Delete(i);
      mObjectList->Insert(name, i-1);
    }

    iter++;
  }
  ObjectSelected(event);
}


void FrameMain::MoveObjecsDown( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    mGLViewer->MoveObject((*iter)->GetName(), 1);

    wxString name((*iter)->GetName().c_str(), wxConvUTF8);
    int i = mObjectList->FindString(name);
    if (i < mObjectList->GetCount()-1) {
      mObjectList->Delete(i);
      mObjectList->Insert(name, i+1);
    }

    iter++;
  }
  ObjectSelected(event);
}


void FrameMain::TextCtrlFocus( wxFocusEvent& event )
{
  wxTextCtrl * ctrl = dynamic_cast<wxTextCtrl *>(event.GetEventObject());
  ctrl->SetSelection(-1,-1);
}


void FrameMain::TransformObjects(wxCommandEvent& event)
{
  double scaleX, scaleY, scaleZ;
  double translateX, translateY, translateZ;
  double rotateX, rotateY, rotateZ;

  if (!mScaleX->GetValue().ToDouble(&scaleX))  scaleX = 1;
  if (!mScaleY->GetValue().ToDouble(&scaleY))  scaleY = 1;
  if (!mScaleZ->GetValue().ToDouble(&scaleZ))  scaleZ = 1;

  if (!mTranslateX->GetValue().ToDouble(&translateX))  translateX = 0;
  if (!mTranslateY->GetValue().ToDouble(&translateY))  translateY = 0;
  if (!mTranslateZ->GetValue().ToDouble(&translateZ))  translateZ = 0;

  if (!mRotateX->GetValue().ToDouble(&rotateX))  rotateX = 0;
  if (!mRotateY->GetValue().ToDouble(&rotateY))  rotateY = 0;
  if (!mRotateZ->GetValue().ToDouble(&rotateZ))  rotateZ = 0;

  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Geometry * object = dynamic_cast<Geometry *>(*iter);
    if (object == NULL) {
      std::cerr << "Warning: Object '" << (*iter)->GetName() << "' is not a geometry - can't be transformed" << std::endl;
    } else {
      object->Scale(scaleX, scaleY, scaleZ);
      object->Translate(translateX, translateY, translateZ);
      object->Rotate(rotateX*M_PI/180.0, rotateY*M_PI/180.0, rotateY*M_PI/180.0);
      UpdateDependentObjects(object);
    }
    iter++;
  }

  mGLViewer->Render();
}


void FrameMain::VisualizeWireframe( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    (*iter)->SetWireframe(event.IsChecked());
    iter++;
  }

  mGLViewer->Render();
}


void FrameMain::VisualizeMeshNormals( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    (*iter)->SetShowNormals(event.IsChecked());
    iter++;
  }

  mGLViewer->Render();
}


void FrameMain::OpacityChanged( wxScrollEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    (*iter)->SetOpacity(event.GetInt()/100.0f);
    iter++;
  }

  mGLViewer->Render();
}


void FrameMain::SetColormap( wxCommandEvent& event )
{
  ColorMap * map = ColorMapFactory::New(std::string(event.GetString().mb_str()));
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    (*iter)->mAutoMinMax = mAutoMinMax->IsChecked();
    double tmp1 = 0;
    mMin->GetValue().ToDouble(&tmp1);
    (*iter)->mMinCMap = tmp1;
    double tmp2 = 1;
    mMax->GetValue().ToDouble(&tmp2);
    (*iter)->mMaxCMap = tmp2;

    (*iter)->SetColorMap(map);
    iter++;
  }

  mColorMapChoice->SetSelection(0);
  mGLViewer->Render();

  /*
  if(!objects.empty()){
    wxString min, max;
    min.Printf(_T("%.3f"), objects.front()->mMinCMap);
    max.Printf(_T("%.3f"), objects.front()->mMaxCMap);
    mMin->SetValue(min);
    mMax->SetValue(max);
    mGLViewer->Render();
  }
  */
}


void FrameMain::SetVisualizationMode( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  int selection = mVisualizationModeChoice->GetSelection();
  if (selection == 0) return;

  VisualizationModeData * mode = dynamic_cast<VisualizationModeData *>(mVisualizationModeChoice->GetClientObject(selection));
  while (iter != iend) {
    (*iter)->SetVisualizationMode(mode->GetVisualizationMode());
    iter++;
  }

  mVisualizationModeChoice->SetSelection(0);
  mGLViewer->Render();
}


void FrameMain::ScaleChanged( wxCommandEvent& event )
{
  if (mUniformScaling->IsChecked()) {
    mScaleY->SetValue(mScaleX->GetValue());
    mScaleZ->SetValue(mScaleX->GetValue());
  }
}


void FrameMain::ToggleUniformScaling( wxCommandEvent& event )
{
  mScaleY->Enable(!mUniformScaling->IsChecked());
  mScaleZ->Enable(!mUniformScaling->IsChecked());
}


void FrameMain::ToggleAutoMinMax( wxCommandEvent& event )
{
  mMin->Enable(!mAutoMinMax->IsChecked());
  mMax->Enable(!mAutoMinMax->IsChecked());
}


void FrameMain::SaveMesh( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {

    // Try to fetch a mesh...
    Mesh * mesh = dynamic_cast<Mesh *>(*iter);
    //Lab4->
    if (mesh == NULL) {
      Implicit * impl = dynamic_cast<Implicit *>(*iter);
      if (impl != NULL)
        mesh = impl->GetMesh();
    }
    //Lab4<-

    if (mesh != NULL) {
      wxFileDialog * dialog = new wxFileDialog(this,_T("Save mesh '") + wxString(mesh->GetName().c_str(), wxConvUTF8) + _T("' as"),_T("."),wxString(mesh->GetName().c_str(), wxConvUTF8) + _T(".obj"),_T("OBJ (*.obj)|*.obj"),wxFD_SAVE, wxDefaultPosition);
      if (dialog->ShowModal() == wxID_OK) {
        wxString filename = dialog->GetPath();
        std::ofstream out(filename.mb_str());
        mesh->save(out);
      }
    }
    iter++;
  }
}


void FrameMain::CaptureScreen( wxCommandEvent& event )
{
	//Bug in windows. ShowModal doesn't work.
	/*
   wxTextEntryDialog filenameDialog(this, _T("Filename"), _T("Enter filename"), _T("dump.tga"));
   wxString filename;
   if (filenameDialog.ShowModal() == wxID_OK) {
   filename = filenameDialog.GetValue();
   wxTextEntryDialog magnDialog(this, _T("Magnification"), _T("Enter magnification"), _T("1"));
   wxString magn;
   if (magnDialog.ShowModal() == wxID_OK) {
   magn = magnDialog.GetValue();
   }

   if (filename == _T("") || magn == _T(""))
   return;

   double mag;
   if (!magn.ToDouble(&mag)) return;
   mGLViewer->ScreenCapture(std::string(filename.mb_str()), mag);
   }*/
	wxFileDialog * dialog = new wxFileDialog(this,_T("Save as"),_T("."),_T("dump.tga"),_T("TGA (*.tga)|*.tga"),wxFD_SAVE, wxDefaultPosition);
	if (dialog->ShowModal() == wxID_OK) {
    wxString filename = dialog->GetPath();
    mGLViewer->ScreenCapture(std::string(filename.mb_str()), 1.0);
	}

}


void FrameMain::Dilate( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Geometry * object = dynamic_cast<Geometry *>(*iter);
    if (object == NULL) {
      std::cerr << "Warning: Object '" << (*iter)->GetName() << "' is not a geometry - can't be transformed" << std::endl;
    } else {
      object->Dilate(GetAmount());
      UpdateDependentObjects(object);
    }
    iter++;
  }

  mGLViewer->Render();
}


void FrameMain::Erode( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Geometry * object = dynamic_cast<Geometry *>(*iter);
    if (object == NULL) {
      std::cerr << "Warning: Object '" << (*iter)->GetName() << "' is not a geometry - can't be transformed" << std::endl;
    } else {
      object->Erode(GetAmount());
      UpdateDependentObjects(object);
    }
    iter++;
  }

  mGLViewer->Render();
}


void FrameMain::Smooth( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Geometry * object = dynamic_cast<Geometry *>(*iter);
    if (object == NULL) {
      std::cerr << "Warning: Object '" << (*iter)->GetName() << "' is not a geometry - can't be transformed" << std::endl;
    } else {
      object->Smooth(GetAmount());
      UpdateDependentObjects(object);
    }
    iter++;
  }

  mGLViewer->Render();
}


double FrameMain::GetAmount()
{
  double amount;
  if (!mAmount->GetValue().ToDouble(&amount)) {
    amount = 1;
    mAmount->SetValue(_T("1"));
  }

  return amount;
}


void FrameMain::HideAllPanels()
{
  mPanelTransform->Hide();
  mPanelVisualization->Hide();
  mPanelDecimation->Hide();
  mPanelSubdivision->Hide();
  mPanelImplicit->Hide();
  mPanelLevelset->Hide();
  mPanelFluid->Hide();

  mPanelSideBar->FitInside();
}


void FrameMain::UpdatePanels()
{
  HideAllPanels();

  // Clear the visualization mode choice box
  while (mVisualizationModeChoice->GetCount() > 1)
    mVisualizationModeChoice->Delete(mVisualizationModeChoice->GetCount()-1);

  // Get selected objects
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {

    // Activate all panels connected to this class
    std::list<wxWindow *>::iterator panelIter = mPanelSwitches[(*iter)->GetTypeName()].begin();
    std::list<wxWindow *>::iterator panelIend = mPanelSwitches[(*iter)->GetTypeName()].end();
    while (panelIter != panelIend) {
      (*panelIter)->Show();
      panelIter++;
    }

    // Add all the visualization modes available for this class
    std::list<GLObject::VisualizationMode> modes = (*iter)->GetVisualizationModes();
    std::list<GLObject::VisualizationMode>::iterator modeIter = modes.begin();
    std::list<GLObject::VisualizationMode>::iterator modeIend = modes.end();
    while (modeIter != modeIend) {
      mVisualizationModeChoice->Append(wxString((*modeIter).GetName().c_str(), wxConvUTF8), new VisualizationModeData(*modeIter));
      modeIter++;
    }

    mVisualizeOpacity->SetValue((*iter)->GetOpacity()*100);
    iter++;
  }

  mPanelSideBar->FitInside();
}





//Lab1->

void FrameMain::AddObjectSimpleMesh( wxCommandEvent& event )
{
  wxFileDialog * dialog = new wxFileDialog(this);
  if (dialog->ShowModal() == wxID_OK) {
    wxString path = dialog->GetPath();
    AddMesh<SimpleMesh>(path);
  }

  delete dialog;
  mGLViewer->Render();
}


void FrameMain::AddObjectHalfEdgeMesh( wxCommandEvent& event )
{
  wxFileDialog * dialog = new wxFileDialog(this);
  if (dialog->ShowModal() == wxID_OK) {
    wxString path = dialog->GetPath();
    AddMesh<HalfEdgeMesh>(path);
  }
  delete dialog;
  mGLViewer->Render();
}


//Lab1<-

//Lab2->

void FrameMain::AddObjectSimpleDecimationMesh( wxCommandEvent& event )
{
  wxFileDialog * dialog = new wxFileDialog(this);
  if (dialog->ShowModal() == wxID_OK) {
    wxString path = dialog->GetPath();
    AddMesh<SimpleDecimationMesh>(path);
  }
  delete dialog;
  mGLViewer->Render();
}


void FrameMain::AddObjectQuadricDecimationMesh( wxCommandEvent& event )
{
  wxFileDialog * dialog = new wxFileDialog(this);
  if (dialog->ShowModal() == wxID_OK) {
    wxString path = dialog->GetPath();
    AddMesh<QuadricDecimationMesh>(path);
  }
  delete dialog;
  mGLViewer->Render();
}


void FrameMain::DecimateObjects( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    DecimationMesh * mesh = dynamic_cast<DecimationMesh *>(*iter);
    if (mesh == NULL)
      std::cerr << "Warning: Object '" << (*iter)->GetName() << "' is not a decimation mesh - can't be decimated" << std::endl;
    else
    {
      /*mesh->decimate();*/
      long targetFaces;
      if (m_DecimationTargetTxtBox->GetValue().ToLong(&targetFaces)) {
        std::cout << "Target Decimation Faces: " << targetFaces << std::endl;
        mesh->decimate( targetFaces );
      }
      else
        std::cout << "Decimating one edge" << std::endl;
      mesh->decimate();
    }

    iter++;
  }

  mGLViewer->Render();
}


//Lab2<-

//Lab3->

void FrameMain::AddObjectCubicSpline( wxCommandEvent& event )
{
  std::vector<Vector3<float> > c;
  c.push_back(Vector3<float>(-.8, 1, .5));
  c.push_back(Vector3<float>(-.7, .8, .5));
  c.push_back(Vector3<float>(-.5, .7, .5));
  c.push_back(Vector3<float>(-.5, 1, .5));
  c.push_back(Vector3<float>(0, .1, .5));
  c.push_back(Vector3<float>(.3, .6, .5));
  c.push_back(Vector3<float>(.4, .5, .5));
  c.push_back(Vector3<float>(.5, 1, .5));

  UniformCubicSpline * curve = new UniformCubicSpline(c, Vector3<float>(1, 0, 0));
  curve->SetName("Cubic spline");
  AddUniqueObject(curve);
  mGLViewer->Render();
}


void FrameMain::AddObjectSubdivisionCurve( wxCommandEvent& event )
{
  std::vector<Vector3<float> > c;
  c.push_back(Vector3<float>(-.8, 1, .5));
  c.push_back(Vector3<float>(-.7, .8, .5));
  c.push_back(Vector3<float>(-.5, .7, .5));
  c.push_back(Vector3<float>(-.5, 1, .5));
  c.push_back(Vector3<float>(0, .1, .5));
  c.push_back(Vector3<float>(.3, .6, .5));
  c.push_back(Vector3<float>(.4, .5, .5));
  c.push_back(Vector3<float>(.5, 1, .5));

  UniformCubicSplineSubdivisionCurve * curve = new UniformCubicSplineSubdivisionCurve(c);
  curve->SetName("Subdivision curve");
  AddUniqueObject(curve);
  mGLViewer->Render();
}


void FrameMain::AddObjectLoopSubdivisionMesh( wxCommandEvent& event )
{
  wxFileDialog * dialog = new wxFileDialog(this);
  if (dialog->ShowModal() == wxID_OK) {
    wxString path = dialog->GetPath();
    AddMesh<LoopSubdivisionMesh>(path);
  }
  delete dialog;
  mGLViewer->Render();
}


void FrameMain::AddObjectStrangeSubdivisionMesh( wxCommandEvent& event )
{
  wxFileDialog * dialog = new wxFileDialog(this);
  if (dialog->ShowModal() == wxID_OK) {
    wxString path = dialog->GetPath();
    AddMesh<StrangeSubdivisionMesh>(path);
  }
  delete dialog;
  mGLViewer->Render();
}


void FrameMain::SubdivideObjects( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Subdivision * subdiv = dynamic_cast<Subdivision *>(*iter);
    if (subdiv == NULL)
      std::cerr << "Warning: Object '" << (*iter)->GetName() << "' is not a subdivision object - can't be subdivided" << std::endl;
    else
    {
      std::cout << "Subdividing" << std::endl;
      subdiv->Subdivide();
    }

    iter++;
  }

  mGLViewer->Render();
}


//Lab3<-

//Lab4->

void FrameMain::AddObjectImplicitSphere( wxCommandEvent& event )
{
  Sphere * sphere = new Sphere(1);
  sphere->SetMeshSampling(GetMeshSampling());
  sphere->SetDifferentialScale(GetDifferentialScale());
  sphere->Triangulate<SimpleMesh>();

  std::cout << "Volume of sphere: " << sphere->ComputeVolume() << std::endl;

  sphere->SetName("Sphere");
  AddUniqueObject(sphere);
  mGLViewer->Render();
}


void FrameMain::AddObjectImplicitMesh( wxCommandEvent& event )
{
  wxFileDialog * dialog = new wxFileDialog(this);
  if (dialog->ShowModal() == wxID_OK) {
    wxString path = dialog->GetPath();
    wxString filename = path.AfterLast('/');
    if (filename == path) // If we're on Windows
      filename = path.AfterLast('\\');
    wxString suffix = path.AfterLast('.');

    if (suffix == _T("obj")) {
      // Create new mesh
      SimpleMesh * mesh = new SimpleMesh();

      // Load mesh
      std::ifstream infile;
      ObjIO objIO;
      infile.open(path.mb_str());
      objIO.Load(mesh, infile);

      // Create new implicit mesh with loaded mesh as argument
      ImplicitMesh * implicitMesh = new ImplicitMesh(mesh);
      implicitMesh->SetDifferentialScale(GetDifferentialScale());
      implicitMesh->SetMeshSampling(GetMeshSampling());
      implicitMesh->Triangulate<SimpleMesh>();

      implicitMesh->SetName("Implicit " + std::string(filename.mb_str()));
      AddUniqueObject(implicitMesh);
    }
    else
      std::cerr << "Error: File type not supported" << std::endl;
  }

  delete dialog;
  mGLViewer->Render();
}


void FrameMain::AddObjectQuadricPlane( wxCommandEvent& event )
{
  Matrix4x4<float> M;
  // Construct the quadric matrix here

  Quadric * Q = new Quadric(M);
  Q->SetBoundingBox(Bbox(-1,1));
  Q->SetMeshSampling(GetMeshSampling());
  Q->SetDifferentialScale(GetDifferentialScale());
  Q->Triangulate<SimpleMesh>();
  Q->SetName("Plane");
  AddUniqueObject(Q);
  mGLViewer->Render();
}


void FrameMain::AddObjectQuadricCylinder( wxCommandEvent& event )
{
  Matrix4x4<float> M;
  // Construct the quadric matrix here

  Quadric * Q = new Quadric(M);
  Q->SetBoundingBox(Bbox(-1,1));
  Q->SetMeshSampling(GetMeshSampling());
  Q->SetDifferentialScale(GetDifferentialScale());
  Q->Triangulate<SimpleMesh>();
  Q->SetName("Cylinder");
  AddUniqueObject(Q);
  mGLViewer->Render();
}


void FrameMain::AddObjectQuadricEllipsoid( wxCommandEvent& event )
{
  Matrix4x4<float> M;
  // Construct the quadric matrix here

  Quadric * Q = new Quadric(M);
  Q->SetBoundingBox(Bbox(-1,1));
  Q->SetMeshSampling(GetMeshSampling());
  Q->SetDifferentialScale(GetDifferentialScale());
  Q->Triangulate<SimpleMesh>();
  Q->SetName("Ellipsoid");
  AddUniqueObject(Q);
  mGLViewer->Render();
}


void FrameMain::AddObjectQuadricCone( wxCommandEvent& event )
{
  Matrix4x4<float> M;
  // Construct the quadric matrix here

  Quadric * Q = new Quadric(M);
  Q->SetBoundingBox(Bbox(-1,1));
  Q->SetMeshSampling(GetMeshSampling());
  Q->SetDifferentialScale(GetDifferentialScale());
  Q->Triangulate<SimpleMesh>();
  Q->SetName("Cone");
  AddUniqueObject(Q);
  mGLViewer->Render();
}


void FrameMain::AddObjectQuadricParaboloid( wxCommandEvent& event )
{
  Matrix4x4<float> M;
  // Construct the quadric matrix here

  Quadric * Q = new Quadric(M);
  Q->SetBoundingBox(Bbox(-2,2));
  Q->SetMeshSampling(GetMeshSampling());
  Q->SetDifferentialScale(GetDifferentialScale());
  Q->Triangulate<SimpleMesh>();
  Q->SetName("Paraboloid");
  AddUniqueObject(Q);
  mGLViewer->Render();
}


void FrameMain::AddObjectQuadricHyperboloid( wxCommandEvent& event )
{
  Matrix4x4<float> M;
  // Construct the quadric matrix here

  Quadric * Q = new Quadric(M);
  Q->SetBoundingBox(Bbox(-2,2));
  Q->SetMeshSampling(GetMeshSampling());
  Q->SetDifferentialScale(GetDifferentialScale());
  Q->Triangulate<SimpleMesh>();
  Q->SetName("Hyperboloid");
  AddUniqueObject(Q);
  mGLViewer->Render();
}


void FrameMain::AddObjectScalarCutPlane( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Implicit * impl = dynamic_cast<Implicit *>(*iter);
    if (impl != NULL) {
      ImplicitValueField * field = new ImplicitValueField(impl);
      ScalarCutPlane * plane = new ScalarCutPlane("Scalar cut plane", 0.005, field);
      AddUniqueObject(plane);

      // Scale cut plane to fill the bounding box
      const Bbox & b = impl->GetBoundingBox();
      float x = b.pMax[0] - b.pMin[0];
      float y = b.pMax[1] - b.pMin[1];
      float z = b.pMax[2] - b.pMin[2];
      float scale = x;
      if (scale < y)  scale = y;
      if (scale < z)  scale = z;
      plane->Scale(scale*0.5);

      mDependentObjects[impl->GetName()].push_back(plane->GetName());
    }

    iter++;
  }
  mGLViewer->Render();
}


void FrameMain::AddObjectVectorCutPlane( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Implicit * impl = dynamic_cast<Implicit *>(*iter);
    if (impl != NULL) {
      ImplicitGradientField * field = new ImplicitGradientField(impl);
      VectorCutPlane * plane = new VectorCutPlane("Vector cut plane", 0.1, field);
      AddUniqueObject(plane);

      // Scale cut plane to fill the bounding box
      const Bbox & b = impl->GetBoundingBox();
      float x = b.pMax[0] - b.pMin[0];
      float y = b.pMax[1] - b.pMin[1];
      float z = b.pMax[2] - b.pMin[2];
      float scale = x;
      if (scale < y)  scale = y;
      if (scale < z)  scale = z;
      plane->Scale(scale*0.5);

      mDependentObjects[impl->GetName()].push_back(plane->GetName());
    }

    iter++;
  }
  mGLViewer->Render();
}


void FrameMain::Union( wxCommandEvent& event )
{
  Implicit * impl = CSG< ::Union, ::BlendedUnion>("U");
  impl->SetMeshSampling(GetMeshSampling());
  impl->SetDifferentialScale(GetDifferentialScale());
  impl->Triangulate<SimpleMesh>();
  AddUniqueObject(impl);
  mGLViewer->SelectObject(impl->GetName());
  mGLViewer->Render();
}


void FrameMain::Intersection( wxCommandEvent& event )
{
  Implicit * impl = CSG< ::Intersection, ::BlendedIntersection>("A");
  impl->SetMeshSampling(GetMeshSampling());
  impl->SetDifferentialScale(GetDifferentialScale());
  impl->Triangulate<SimpleMesh>();
  AddUniqueObject(impl);
  mGLViewer->SelectObject(impl->GetName());
  mGLViewer->Render();
}


void FrameMain::Difference( wxCommandEvent& event )
{
  Implicit * impl = CSG< ::Difference, ::BlendedDifference>("-");
  impl->SetMeshSampling(GetMeshSampling());
  impl->SetDifferentialScale(GetDifferentialScale());
  impl->Triangulate<SimpleMesh>();
  AddUniqueObject(impl);
  mGLViewer->SelectObject(impl->GetName());
  mGLViewer->Render();
}


void FrameMain::SwitchBlending( wxCommandEvent& event )
{
  mBlendLabel->Enable(event.IsChecked());
  mBlendParameter->Enable(event.IsChecked());
}


void FrameMain::ResampleImplicit( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Implicit * impl = dynamic_cast<Implicit *>(*iter);
    if (impl != NULL) {
      impl->SetMeshSampling(GetMeshSampling());
      impl->Update();
      impl->Triangulate<SimpleMesh>();
    }
    iter++;
  }
  mGLViewer->Render();
}


void FrameMain::DifferentialScaleChanged( wxScrollEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();
  std::list<GLObject *>::iterator iter = objects.begin();
  std::list<GLObject *>::iterator iend = objects.end();
  while (iter != iend) {
    Implicit * impl = dynamic_cast<Implicit *>(*iter);
    if (impl != NULL) {
      impl->SetDifferentialScale(GetDifferentialScale());
      impl->Update();
    }
    iter++;
  }
  mGLViewer->Render();
}

double FrameMain::GetMeshSampling()
{
  double sampling;
  if (!mMeshSampling->GetValue().ToDouble(&sampling)) {
    sampling = 0.05;
    mMeshSampling->SetValue(_T("0.05"));
  }

  return sampling;
}

double FrameMain::GetDifferentialScale()
{
  int scale = mDifferentialScale->GetValue();
  if (scale == 0)
    scale = 1;
  return scale/100.0f;
}


//Lab4<-


