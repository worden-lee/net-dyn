#include "DotDisplay.h"

DotGraphController::DotGraphController(double dperiod, double fperiod)
  : DisplayController<DotDisplay>(dperiod, fperiod)
{}

void DotGraphController::createDisplay()
{ display = new DotDisplay(this); }

void DotGraphController::_update(void)
{
  if (allowRecording())
    record();
  if (allowDisplaying())
    DisplayController<DotDisplay>::updateDisplay();
}
void DotGraphController::record()
{
  if (allowRecording())
  { if (!display)
      createDisplay();
    display->record();
  }
}

void DotGraphController::updateDisplay(void)
{
  _update();
}

DotGraphController::~DotGraphController()
{ delete display; }
