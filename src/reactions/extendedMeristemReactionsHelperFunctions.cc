﻿
// Filename     : extendedMeristemReactionsHelperFunctions.cc
// Description  : Helper Functions for extended Meristem sims
// Author(s)    : Al Do (ado012@ucr.edu)
// Created      : June 2019
// Revision     : $Id: grn.cc 669 2016-05-12 08:53:12Z korsbo $
//

#include"baseReaction.h"
#include"extendedMeristemReactionsHelperFunctions.h"
#include"../organism.h"
#include <random> //ADDITION 051517 Currently unusable due to problems with C11
#include <fstream>//ADDITION 051517
#include <stdlib.h>//ADDITION 051517
#include <string> //ADDITION 100218






double Clavata3ActivationMechanisms(int clv3ActivateMethod, double clv3Creation, double clv3P, int wusMonomer,
                                    int wusDimer, double wusMonomerCoefficient, double wusDimerCoefficient,
                                    int chromoCycle, Compartment &compartment, double polTimeLimit, int crmOrMarkerSwitch,
                                    double distanceFromCentralAxis, double clv3SourceWidth, double zvar, double clv3ZAxisLimit) //Change: added new parameters to bring distance finding inside 122821
{

    //Clavata3ActivationMechanisms: Determines Activation based upon a selected strategy. Currently using method 9



    //Activation strategies
    if (clv3ActivateMethod==1)//graded monomer activation/graded dimer repression
    {
        //clv3Creation=clv3P*(wusMonomerCoefficient*wusMonomer-wusDimerCoefficient*wusDimer);
        //alter activation function to make the crm config and not just the wus concentration generated by config play a role in determining clv3 output
        clv3Creation=clv3P*(wusMonomerCoefficient*wusMonomer-wusDimerCoefficient*wusDimer);


        if (clv3Creation < 0)
            clv3Creation=0;
    }
    else if (clv3ActivateMethod==2)//single monomer activation
    {
        if (wusMonomer>=1)
            clv3Creation=clv3P*5*(wusMonomerCoefficient*1);
    }
    else if (clv3ActivateMethod==3)//graded monomer activation/ single dimer repression-dimer dominant
    {
        clv3Creation=clv3P*(wusMonomerCoefficient*wusMonomer-wusDimerCoefficient*wusDimer);
        if (clv3Creation < 0)
            clv3Creation=0;
    }
    else if (clv3ActivateMethod==4) //graded monomer activation
        clv3Creation=clv3P*(wusMonomerCoefficient*wusMonomer);
    else if (clv3ActivateMethod==5) //monomer/(dimer+1)
        clv3Creation=clv3P*(wusMonomerCoefficient*wusMonomer/((wusDimerCoefficient*wusDimer)+1));

    else if (clv3ActivateMethod==6)//monomer more than half
    {
        if (wusMonomer>=3)
            clv3Creation=clv3P*(5*wusMonomerCoefficient*wusMonomer);
    }
    else if (clv3ActivateMethod==7)//single monomer single dimer -dimer dominant
    {
        if (wusMonomer>=1 && wusDimer==0)
            clv3Creation=clv3P*(5*wusMonomerCoefficient*wusMonomer);
    }
    else if (clv3ActivateMethod==8)//polymerase activation
    {


        if (compartment.polymeraseBindSite[crmOrMarkerSwitch][chromoCycle]==1 && (compartment.rnaPolTimer[crmOrMarkerSwitch][chromoCycle] <= 0))
       {
            //will have to include this in other mechanisms if you use them in the future 122821
            //only create clv3 in the clv3 generation range
            if (distanceFromCentralAxis<=clv3SourceWidth && zvar>=clv3ZAxisLimit)
            {
            clv3Creation=1*clv3P;

            //tracking statement
            compartment.polFire +=1;
            }

            else
            {
             clv3Creation=0;
            }

        //unbind polymerase and zero timer
        compartment.polymeraseBindSite[crmOrMarkerSwitch][chromoCycle]=0;
        compartment.rnaPolTimer[crmOrMarkerSwitch][chromoCycle]=0;


        }
    }



    return clv3Creation;
}

//Restrict to internal cylinder through multiplication




void CRMProbabilityGenerator(int m4Flag, int crmOrMarkerSwitch, int chromoCycle, probabilitySegment *probabilityMatrix, int *crmOccupancy, int i, double crmActivityCoefficient,
                             double cooptMonEffect, double cooptDimEffect, double geneCRMSiteBindMaxBaseChance, double geneCRMSiteChance_Unbind, double concModifier, int& eventNum,
                             double& probabilityDeltaSum, int HACisEleBonusCoopt, int neighborOnlyCoopt, double dimerBonusBindP, double polBaseBindAffinity, Compartment &compartment, double dimerUnbindP,
                             int L1nodimer, int bonusL1MonCoopt, double distanceFromBase)
{

    //1: CRMProbabilityGenerator: For each site calculate event probabilities based on binding constant, site state, concentration, and site's base affinity
    //build a probability matrix to allow for an event to be chosen at each timestep.




    //2: Setup
    //0->1
    //1->0 1->2
    //2->1
    int neighbors=0;
    int nearNeighbor1 =0;
    int nearNeighbor2 = 0;
    int nearNeighbor1Location = 0;
    int nearNeighbor2Location = 0;
    int neighborsInInterval = 0;
    int cooperative970Bonus=0;
    double cooperativeEffect, eventChance, eventChance2, polBindChance;
    //int *arraysize=crmOccupancy.size();
    int sampArray[5]={crmOccupancy[0],crmOccupancy[1],crmOccupancy[2],crmOccupancy[3],crmOccupancy[4]};
    int chromoCycleBonus=0;

    if (crmOrMarkerSwitch == 1)
        chromoCycleBonus = 2;


//tag the probability segment with the chromosome it is referring to.
    //probabilityMatrix[eventNum].chromosome=chromoCycle; source of probabilityMatrix tagging bug?



    /*
    //ADDITION 070322: AD Scale unbind to bind rates other than 0.1
       double unbindmultiplier = 0;

       unbindmultiplier = geneCRMSiteBindMaxBaseChance/0.1;

       geneCRMSiteChance_Unbind = geneCRMSiteChance_Unbind * unbindmultiplier;
       //ADDITION 070322
*/



    //3: Empty Site Probability Generation: if site is empty, a binding event can occur
    if (crmOccupancy[i]==0)
    {//3b: calculate chance to bind
        eventChance=crmActivityCoefficient*geneCRMSiteBindMaxBaseChance*concModifier;


        if (eventChance > 0)//3c: set event on probability spectrum
        {
            probabilityMatrix[eventNum].begin=probabilityDeltaSum;
            probabilityDeltaSum +=eventChance;
            probabilityMatrix[eventNum].end=probabilityDeltaSum;
            probabilityMatrix[eventNum].site=i;
            probabilityMatrix[eventNum].action="BIND";
            probabilityMatrix[eventNum].chromosome=chromoCycle;
            eventNum++;


        }


    }



    //4: Monomer Occupied Probability Generation: if site is monomer occupied, a binding or unbinding event can occur.
    else if (crmOccupancy[i]==1 )
    {//4b: Bind Calculation: calculate chance to bind
        eventChance= crmActivityCoefficient*geneCRMSiteBindMaxBaseChance*concModifier*dimerBonusBindP;
//if L1nodimer is activated. No dimerization in L1
        if ((eventChance > 0 && L1nodimer == 0) || (eventChance > 0 && L1nodimer == 1 && distanceFromBase < 8.5))
        {//4c: add bind event to probability matrix
            probabilityMatrix[eventNum].begin=probabilityDeltaSum;
            probabilityDeltaSum +=eventChance;
            probabilityMatrix[eventNum].end=probabilityDeltaSum;
            probabilityMatrix[eventNum].site=i;
            probabilityMatrix[eventNum].action="BIND";
            probabilityMatrix[eventNum].chromosome=chromoCycle;

            eventNum++;
        }


        //4d: Neighbor calculation: neighboring monomers reduce chance of monomer unbinding. Currently using option 2

        if (neighborOnlyCoopt==1)
        {
            for (int m=0; m<=4; m++)
            {if (crmOccupancy[m]==1 && m !=i && (m==(i-1) || m==(i+1)))
                 neighbors += 1; }

            if (crmOccupancy[1]==1 and m4Flag==1 && (i==0 || i==2)) //if 970 has monomer, give cooperative bonus
                cooperative970Bonus=HACisEleBonusCoopt;



        }

        else if (neighborOnlyCoopt==2)//Betul style cooperativity: Monomer
        {

            if (crmOccupancy[1]==1) //if 970 occupied by monomer, cooperativity is possible
            {

                if (i == 1) //if examined site is 970, look for neighbors
                {
                    for (int m=0; m<=4; m++)
                    {

                        if (crmOccupancy[m]==1 && m != 1)//if site is monomer occupied and not 970 add it as a neighbor
                            neighbors += 1;

                    }


                }


                else  //if you are looking at another site and 970 is occupied. It counts as 1 neighbor
                {
                neighbors += 1;

                }


            }


        }


        else
        {
            for (int m=0; m<=4; m++)
            {if (crmOccupancy[m]==1 && m !=i)
                 neighbors += 1; }
        }




        //4e: Bonus Cooperativity calculation
        if (crmOccupancy[1]==1) //bonus cooperativity if 970 is occupied by monomer, not just 970m4 according to Betul.can it get bonus even if its by itself?
        cooperative970Bonus=HACisEleBonusCoopt;


        //NOT USED: If bonus L1 Monomer Cooperativity is active add cooperativity
        if (crmOccupancy[1]==1 && bonusL1MonCoopt != 0)
            cooperative970Bonus += bonusL1MonCoopt;


        //4f: Total Cooperativity calculation
        //added monomer only coopt 122020
        cooperativeEffect = pow(cooptMonEffect,(neighbors+cooperative970Bonus));
        cooperative970Bonus=0;


        //4g: Unbind Probability Calculation: calculate chance to unbind (check that cooperative effect doesn't kick in if neighbors = 0)
        if (neighbors>0)
        {
            eventChance2=crmActivityCoefficient*geneCRMSiteChance_Unbind*cooperativeEffect;

        }

        else
            eventChance2=crmActivityCoefficient*geneCRMSiteChance_Unbind;

        neighbors=0;
        cooperativeEffect=0; //zero out cooperativeEffect to prevent carry over from other sites


        //4f: Add Unbind Probability to matrix
        if (eventChance2 > 0)
        {
            probabilityMatrix[eventNum].begin=probabilityDeltaSum;
            probabilityDeltaSum += eventChance2;
            probabilityMatrix[eventNum].end=probabilityDeltaSum;
            probabilityMatrix[eventNum].site=i;
            probabilityMatrix[eventNum].action="UNBIND";
            probabilityMatrix[eventNum].chromosome=chromoCycle;

            eventNum++;
        }


        //4g: Calculate and add RNA polymerase bind event
        if (compartment.polymeraseBindSite[crmOrMarkerSwitch][chromoCycle]==0)
        {
            //Doesn't need specific monomers because probability is automatically increased the more sites are involved
            polBindChance=crmActivityCoefficient*polBaseBindAffinity;
            probabilityMatrix[eventNum].begin=probabilityDeltaSum;
            probabilityDeltaSum +=polBindChance;
            probabilityMatrix[eventNum].end=probabilityDeltaSum;
            probabilityMatrix[eventNum].site=i;
            probabilityMatrix[eventNum].action="POLBIND";
            probabilityMatrix[eventNum].chromosome=chromoCycle;
            eventNum++;
        }

    }

    //5: Dimer Probability Generation: if site is dimer occupied, an unbinding event can occur.
    else if (crmOccupancy[i]==2)
    {

        //5b: Neighbor Calculation: Currently Option 2 is used
        if (neighborOnlyCoopt==1)
        {
            //cooperative effect, neighboring dimers reduce chance of dimer unbinding.
            for (int m=0; m<=4; m++)
            {if (crmOccupancy[m]==2 && m !=i && (m==i-1 || m==i+1))
                 neighbors += 1; }


            if (crmOccupancy[1]==2 && m4Flag==1 && (i==0 || i==2)) //if 970 has monomer, give cooperative bonus
                cooperative970Bonus=HACisEleBonusCoopt;
        }

        else if (neighborOnlyCoopt==2)//Betul style cooperativity: Dimer
        {

            if (m4Flag==1 && i==1) //m4 gets cooperativity from all other sites.
            {
                for (int m=0; m<=4; m++)
                {
                    if (crmOccupancy[m] == 2 && m != 1)
                        neighbors += 1;
                }

            }


            else //all other sites get cooperativity from nearest neighbor.
            {

                for (int m=0; m<=4; m++)
                {

                    if ((m < i) && (crmOccupancy[m] != 0))
                    {
                        nearNeighbor1 = crmOccupancy[m];


                    }

                    if ((m > i) && (crmOccupancy[m] != 0) && (nearNeighbor2 == 0 ))
                    {
                        nearNeighbor2 = crmOccupancy[m];


                    }


                }

                if (nearNeighbor1 == 2)
                    neighbors += 1;

                if (nearNeighbor2 == 2)
                    neighbors += 1;

            }


        }



        else
        {
            //cooperative effect, neighboring dimers reduce chance of dimer unbinding.
            for (int m=0; m<=4; m++)
            {if (crmOccupancy[m]==2 && m !=i)
                 neighbors += 1; }


        }

        //5c: Bonus Cooperativity Calculation
        if (crmOccupancy[1]==2 && m4Flag==1 ) //if 970m4 has dimer, give cooperative bonus
            cooperative970Bonus=HACisEleBonusCoopt;



        //5d: total cooperativity calculation
        cooperativeEffect = pow(cooptDimEffect,(neighbors+cooperative970Bonus));

        cooperative970Bonus=0;


        //5e: Unbind Calculation

        if (neighbors>0)
            eventChance=crmActivityCoefficient*geneCRMSiteChance_Unbind*cooperativeEffect*dimerUnbindP;


        else
            eventChance=crmActivityCoefficient*geneCRMSiteChance_Unbind*dimerUnbindP;

        neighbors=0;
        cooperativeEffect=0; //zero out cooperativeEffect to prevent carry over from other sites


        if (eventChance > 0)
        {//5f: add to probability matrix
            probabilityMatrix[eventNum].begin=probabilityDeltaSum;
            probabilityDeltaSum += eventChance;
            probabilityMatrix[eventNum].end=probabilityDeltaSum;
            probabilityMatrix[eventNum].site=i;
            probabilityMatrix[eventNum].action="UNBIND";
            probabilityMatrix[eventNum].chromosome=chromoCycle;
            eventNum++;
        }



    }





}







int CRMEventPicker(int site, int eventFlag, std::string action, double eventBegin, double eventEnd,
                   double wusConc, double randvalue2, int eventNum, double crmTimerLength, double &crmSiteActiveTimer,
                   int &crmSite, Compartment &compartment, int chrom, int monRecruitLimit, int crmOrMarkerSwitch, int polTimeLimit, int monResidenceTime)
{


    //EventPicker: Pick event from probability matrix

    //1. Setup
    int monKeepAgeRoll;
    int neighbors = 0;
    int nearNeighbor1 = 0;
    int nearNeighbor2 = 0;
    int neighborsInInterval = 0;
//std::ofstream cellTrackOutputFile;//file testing

    //2: Event Determination: if random value falls on the event location in the matrix it is chosen.
    //what if randvalue falls on margin?
    if (eventBegin < randvalue2 && eventEnd > randvalue2)
    {

        /*
        if (compartment.index() == 89 || compartment.index() == 375 || compartment.index() == 777 ||
               compartment.index() == 101 || compartment.index() == 241 || compartment.index() == 272 ||
               compartment.index() == 87 || compartment.index() == 421 || compartment.index() == 509)
           {

               cellTrackOutputFile.open ("cellTrackOutputFile" , std::fstream::app);

               cellTrackOutputFile << compartment.index() << " " << compartment.timer << " " << compartment.latestRandTimeStep[0] << " "
                                   << compartment.geneCRM[0][0] << " " << compartment.geneCRM[0][1] << " " << compartment.geneCRM[0][2] << " "
                                   << compartment.geneCRM[0][3] << " " << compartment.geneCRM[0][4]  << " "<< compartment.polymeraseBindSite[0][0]<< " "
                                   << action << " " << compartment.polFire << std::endl;

               cellTrackOutputFile.close();
           }

*/






        //2b: A binding event is selected
        if (action=="BIND" && wusConc >=1 )
        {crmSite +=1;


            //initialize timer if empty site is bound by monomer, currently a dimer forming does not reset the timer.
         if (crmSite==1)
         {

                 compartment.monRecruitCountDown[chrom][site]=monRecruitLimit; //How many times a monomer can recruit a polymerase. currently set to 1000

                 //calculate monomer neighbors
                 if (site == 1) //if examined site is 970, look for nearest neighbors
                 {
                     for (int m=0; m<=4; m++)
                     {

                         if (compartment.geneCRM[chrom][m]==1 && m != 1)
                             neighbors += 1;

                     }
                 }



                 else if ( compartment.geneCRM[chrom][1] == 1)
                 {

                         neighbors += 1;

                 }


//timer only for monomers with no cooperativity DISABLED condition for timer CHANGE: 111023
                 //if (neighbors == 0) //if the site is one of the first three at has monomer neighbors. Start a monomer timer
                     compartment.monResidenceTimer[chrom][site]=monResidenceTime; //monResidenceTime set to 100

    
         }


         eventFlag=1;

        }



        else if (action=="UNBIND")//2c: an unbinding event is chosen
        {
            crmSite -=1; //update crm site


            if (crmSite==0)//zero out timers if binding site is empty
            {

                    compartment.monResidenceTimer[chrom][site]=0;
                    compartment.monRecruitCountDown[chrom][site]=0;

           }


            eventFlag=1;


        }

        else if (action=="POLBIND") //2d: A polymerase bind event is chosen
        {

            //bind polymerase site and set polymerase unbind timer and update monomer polymerase recruitment count
            if (compartment.monRecruitCountDown[chrom][site] > 0 && crmOrMarkerSwitch==0)
            {
                compartment.polymeraseBindSite[crmOrMarkerSwitch][chrom]=1;
                compartment.rnaPolTimer[crmOrMarkerSwitch][chrom]=polTimeLimit;
                compartment.monRecruitCountDown[chrom][site] -= 1;
            }



            eventFlag=1;
        }





    }

    //what happens if the random value lands right on a margin?

//3: return flag on whether an event was chosen
    return eventFlag;



}



void CRMSummer(int crmSite, int j, double &crmActiveTimer, int &wusMonomer, int &wusDimer, double timeStepIncrement)
{
    if (crmSite==1)
    {

        wusMonomer +=1;


    }

    else if (crmSite==2)
        wusDimer +=1;

}


int CRMOccupiedTest(int crmSite)
{
    int occupiedFlag=0;

    if (crmSite==1)
    occupiedFlag=1;

    else if (crmSite==2)
    occupiedFlag=1;

    return occupiedFlag;

}






int TimeStepOverflowHandler(double timeStep, int polTimeLimit, double timeStepOverflow, int crmOrMarkerSwitch,
                           double clv3P, double timeStepRemain, Compartment &compartment, OverFlowResults &overFlowResults1)
{
    //TimeStepOverFlowHandler: handles the effects of stochastic timestep overflow
    //returns CLV3Creation, timeStepRemain
    //alters compartment.rnaPolTimer, compartment.polymeraseBindSite, compartment.stochastictimeStepOverflow

    //1. Setup
    double clv3Creation=0;
    int tfI;
    int timerSet[2]= {0,1};




//2. Large Overflow handling > 1
    //If stochastic timestepoverflow is larger than ode timestep then the entire
    //rest of the ODE timestep will be skipped. Overflow will be subtracted by 1
    //Timers will be updated
    if (timeStepOverflow >= timeStep)//Do not need to add new activation because that presumably happened in previous timestep?
    {
        compartment.stochasticStepOverflow[crmOrMarkerSwitch]=compartment.stochasticStepOverflow[crmOrMarkerSwitch]-1;
        timeStepRemain=0;//timestep is zeroed out by large overflow

        compartment.previousTimeStep[crmOrMarkerSwitch]=compartment.latestTimeStep[crmOrMarkerSwitch];
        compartment.latestTimeStep[crmOrMarkerSwitch]=1;
        //no need to update timers here since overflow from previous timesteps are handled in timestep generation




    }

    //3. Small Timestep handling < 1
    //if timestep overflow is smaller, simply subtract the overflow from the time remaining
    else
    {
        timeStepRemain=timeStepRemain-timeStepOverflow;//should not overflow into the next timestep but might be a good idea to check if overflow is actually carried out
        compartment.stochasticStepOverflow[crmOrMarkerSwitch]=0;
        compartment.previousTimeStep[crmOrMarkerSwitch]=compartment.latestTimeStep[crmOrMarkerSwitch];
        compartment.latestTimeStep[crmOrMarkerSwitch]=timeStepOverflow;

    }

//4: Overflow output
    overFlowResults1.clv3=clv3Creation;
    overFlowResults1.timeStepRemain= timeStepRemain;

    //returns CLV3Activation, newTimestepOverflow, newTimestep

    return 0;
}


double TimeStepGenerator(int crmOrMarkerSwitch, double randvalue1, int polTimeLimit, Compartment &compartment)
{
    //TimeStepGenerator: Generates timestep
    //returns timeStepRemain
    //alters compartment.rnaPolTimer, compartment.monResidenceTimer

    //1. Setup
    double timeStepStochastic=0;
    double timeStepIncrement=0;
    double minTimer=0;
    double timers[14]; //one extra space? Corrected mistake that set timers to int instead of double 122621
    int tgI;
       int timerSet[2]= {0,1};



    compartment.probabilityDeltaSum[crmOrMarkerSwitch]=compartment.probabilityDeltaSum[crmOrMarkerSwitch];//scaling factor ask weitao?

    //2. stochastic timestep generation
    timeStepStochastic=(1/((compartment.probabilityDeltaSum[crmOrMarkerSwitch]))*log(1/(randvalue1)));


    //3: Find the minimum timestep by gathering all unbind times and timestep in order to determine next event to increment towards

    //3b: Gather stochastic timestep and polymerase unbind times
    timers[0]=timeStepStochastic;
    timers[1]=compartment.rnaPolTimer[crmOrMarkerSwitch][0];
    timers[2]=compartment.rnaPolTimer[crmOrMarkerSwitch][1];//Examine later. POSSIBLE ISSUE Fixed 121221

    //3c: Gather monomer unbind times
    //Chromosome 1 monomers
    for (tgI = 3; tgI <= 7; tgI++)
    timers[tgI]=compartment.monResidenceTimer[0][tgI-3];

//chromosome 2 monomers
    for (tgI = 8; tgI <= 12; tgI++)
    timers[tgI]=compartment.monResidenceTimer[1][tgI-8];


//4: Determine the smallest time by looping through the gathered times
minTimer=timers[0];//set timer to stochastictimestep initially, does this get rid of the risk of the timer being stuck at 0 even with existing small timers?

    for (tgI=0; tgI <=12; tgI++)
    {
    if ((timers[tgI] < minTimer) && (timers[tgI] > 0)) //make sure to pick smallest not zero timers
        minTimer = timers[tgI];

    }

    //5. Timestep is the minimum time to next event
    timeStepIncrement = minTimer;



      //6. update timers

    //6b: update RNA polymerase timer
   //if the timestep overflows into the next. This should be okay because the rnapoltimer is already updating
      //at beginning of timestep.
    if (compartment.polymeraseBindSite[crmOrMarkerSwitch][0]==1)
        compartment.rnaPolTimer[crmOrMarkerSwitch][0] -= timeStepIncrement;
    if (compartment.polymeraseBindSite[crmOrMarkerSwitch][1]==1)
        compartment.rnaPolTimer[crmOrMarkerSwitch][1] -= timeStepIncrement;


    //6b: update monomer unbinding timers
    for (tgI = 0; tgI<=4; tgI++)
    {
       //chromosome 1
        if (compartment.monResidenceTimer[0][tgI] > 0 && compartment.geneCRM[0][tgI]==1)
        {
            compartment.monResidenceTimer[0][tgI] -= timeStepIncrement;
//unbind expired monomers
            if (compartment.monResidenceTimer[0][tgI] <= 0)
            compartment.geneCRM[0][tgI] = 0;
        }

        //chromosome 2
        if (compartment.monResidenceTimer[1][tgI] > 0 && compartment.geneCRM[1][tgI]==1)
        {
        compartment.monResidenceTimer[1][tgI] -= timeStepIncrement;
//unbind expired monomers
         if (compartment.monResidenceTimer[1][tgI] <= 0)
        compartment.geneCRM[1][tgI] = 0;
        }






    }

    //7: return timestep
return timeStepIncrement;
}






//ResetProbabilityMatrix(probabilityMatrix)

