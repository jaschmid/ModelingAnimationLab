#include "FrameMain.h"

FrameMain::FrameMain( wxWindow* parent ) : BaseFrameMain( parent )
{
  GLGridPlane * plane = new GLGridPlane("Grid");
  plane->SetDimensions(4.0, 4.0);
  AddUniqueObject(plane);

  //GLAxis * axis = new GLAxis("axis");
  //AddUniqueObject(axis);

  // Connect the "object selected" event triggered by GLViewer
  Connect(wxID_ANY, wxEVT_GL_OBJECT_SELECTED, wxCommandEventHandler(FrameMain::ObjectSelected));

  // Add the panels to show for all corresponding types of GLObjects
  mPanelSwitches[typeid(Mesh).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(Mesh).name()].push_back(mPanelVisualization);
  mPanelSwitches[typeid(DecimationMesh).name()].push_back(mPanelTransform);
  mPanelSwitches[typeid(DecimationMesh).name()].push_back(mPanelVisualization);
  mPanelSwitches[typeid(DecimationMesh).name()].push_back(mPanelDecimation);

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

void FrameMain::SaveSelected( wxCommandEvent& event )
{
  std::list<GLObject *> objects = mGLViewer->GetSelectedObjects();

  if (objects.empty() || objects.size() > 1){
    wxMessageDialog * dialog = new wxMessageDialog(this, _T("You need to select one, but only one, object to save."),
                                                   _T("Nothing selected"), wxOK);
    dialog->ShowModal();
    delete dialog;
    return;
  }

  wxFileDialog * dialog = new wxFileDialog(this, _T("Save to file"), _T(""), _T("mesh.obj"), _T("*.obj"), wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
  if (dialog->ShowModal() == wxID_OK) {
    wxString path = dialog->GetPath();
    std::ofstream out( path.mb_str() );

    Mesh * m = dynamic_cast<Mesh *>(objects.front());
   if(m != NULL){
      m->save(out);
    }

  }

  delete dialog;
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

  if(!objects.empty()){
    wxString min, max;
    min.Printf(_T("%.3f"), objects.front()->mMinCMap);
    max.Printf(_T("%.3f"), objects.front()->mMaxCMap);
    mMin->SetValue(min);
    mMax->SetValue(max);
    mGLViewer->Render();
  }
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


void FrameMain::ScaleChanged( wxCommandEvent& event )
{
  if (mUniformScaling->IsChecked()) {
    mScaleY->SetValue(mScaleX->GetValue());
    mScaleZ->SetValue(mScaleX->GetValue());
  }
}

void FrameMain::ToggleAutoMinMax( wxCommandEvent& event )
{
  mMin->Enable(!mAutoMinMax->IsChecked());
  mMax->Enable(!mAutoMinMax->IsChecked());
}

void FrameMain::ToggleUniformScaling( wxCommandEvent& event )
{
  mScaleY->Enable(!mUniformScaling->IsChecked());
  mScaleZ->Enable(!mUniformScaling->IsChecked());
}

void FrameMain::TextCtrlFocus( wxFocusEvent& event )
{
  wxTextCtrl * ctrl = dynamic_cast<wxTextCtrl *>(event.GetEventObject());
  ctrl->SetSelection(-1,-1);
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

/*
 * This event is triggered when an object is selected in the GUI
 */
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

