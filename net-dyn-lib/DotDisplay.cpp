#include "DotDisplay.h"

DotGraphController::DotGraphController(double dperiod, double fperiod)
  : DisplayController<DotDisplay>(dperiod, fperiod)
{}

void DotGraphController::createDisplay()
{ display = new DotDisplay(this); }

DotGraphController::~DotGraphController()
{ delete display; }
