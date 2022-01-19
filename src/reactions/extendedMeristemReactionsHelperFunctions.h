//
// Filename     : extendedMeristemReactionsHelperFunctions.h
// Description  : Helper functions for extended meristem simulation
// Author(s)    : Al Do (ado012@ucr.edu)
// Created      : June 2017
// Revision     :
//
#ifndef EXTENDEDMERISTEMREACTIONSHELPERFUNCTIONS_H
#define EXTENDEDMERISTEMREACTIONSHELPERFUNCTIONS_H

#include <cmath>
#include <math.h>

#include "baseReaction.h"


struct OverFlowResults {
    double clv3=0;
    double timeStepRemain=0;
};



double Clavata3ActivationMechanisms(int activation, double clv3Creation, double clv3P, int wusMonomer, int wusDimer,
                                    double wusMonomerCoefficient, double wusDimerCoefficient, int chromoCycle, Compartment &compartment,
                                    double polTimeLimit, int crmOrMarkerSwitch,double distanceFromCentralAxis, double clv3SourceWidth, double zvar, double clv3Barrier);


void CRMProbabilityGenerator(int m4Flag, int crmOrMarkerSwitch, int chromoCycle, probabilitySegment *probabilityMatrix, int *crmOccupancy, int i, double crmActivityCoefficient, double cooptMonEffect, double cooptDimEffect,
                             double geneCRMSiteBindMaxBaseChance, double geneCRMSiteChance_Unbind, double concModifier, int &eventNum, double &probabilityDeltaSum, int HABonusCoopt, int neighborOnlyCoopt, double dimerBindP,
                             double polBaseBindAffinity, Compartment &compartment, double dimerUnbindP, int L1nodimer, int bonusL1MonCoopt, double distanceFromBase);



int CRMEventPicker(int site, int eventFlag, std::string action, double eventBegin, double eventEnd, double wusConc, double randvalue2, int eventNum, double crmTimerLength, double &crmSiteActiveTimer,
                   int &crmSite, Compartment &compartment, int chrom, int monFireLimit, int crmOrMarkerSwitch, int polTimeLimit, int unbindLimit);


void CRMSummer(int crmSite, int j, double &crmActiveTimer, int &wusMonomer, int &wusDimer, double timeStepIncrement);

int TimeStepOverflowHandler(double timeStep, int polTimeLimit, double timeStepOverflow, int crmOrMarkerSwitch,
                           double clv3P, double timeStepRemain, Compartment &compartment, OverFlowResults &overFlowResults1);


double TimeStepGenerator(int crmOrMarkerSwitch, double randvalue1, int polTimeLimit, Compartment &compartment);

#endif
