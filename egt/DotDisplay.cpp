#include "DotDisplay.h"

DotGraphController::DotGraphController(double dperiod, double fperiod)
  : DisplayController<DotDisplay>(dperiod, fperiod)
{}

void DotGraphController::createDisplay()
{ display = new DotDisplay(this); }

// void DotGraphController::_update(void)
// {
//   if (allowRecording())
//     recordFile();
//   if (allowDisplaying())
//     DisplayController<DotDisplay>::updateDisplay();
// }

void DotGraphController::recordFile()
{
  if (allowRecording())
  { if (!display)
      createDisplay();
    display->recordFile();
  }
}

void DotGraphController::updateDisplay(void)
{
  //  _update();
  if (allowDisplaying())
    DisplayController<DotDisplay>::updateDisplay();
}

DotGraphController::~DotGraphController()
{ delete display; }
