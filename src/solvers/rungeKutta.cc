
// Filename     : rungeKutta.cc
// Description  : Different Runge Kutta solvers
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
// Created      : August 2006
// Revision     : $Id: rungeKutta.cc 672 2016-08-11 11:32:43Z korsbo $
//
#include <cmath>
#include "rungeKutta.h"

////////////////////////////////////////////////////////////////////////
//                                                                                 //
// Dopr853: 8th order Runge-Kutta with adaptive step-length  //
//                                                                                 //
///////////////////////////////////////////////////////////////////////

Dopr853::Dopr853(Organism *O,std::ifstream &IN)
	:BaseSolver(O,IN)
{
  readParameterFile(IN);
}

const double Dopr853::c2  = 0.526001519587677318785587544488e-01;
const double Dopr853::c3  = 0.789002279381515978178381316732e-01;
const double Dopr853::c4  = 0.118350341907227396726757197510e+00;
const double Dopr853::c5  = 0.281649658092772603273242802490e+00;
const double Dopr853::c6  = 0.333333333333333333333333333333e+00;
const double Dopr853::c7  = 0.25e+00;
const double Dopr853::c8  = 0.307692307692307692307692307692e+00;
const double Dopr853::c9  = 0.651282051282051282051282051282e+00;
const double Dopr853::c10 = 0.6e+00;
const double Dopr853::c11 = 0.857142857142857142857142857142e+00;
const double Dopr853::c14 = 0.1e+00;
const double Dopr853::c15 = 0.2e+00;
const double Dopr853::c16 = 0.777777777777777777777777777778e+00;

const double Dopr853::b1 =   5.42937341165687622380535766363e-2;
const double Dopr853::b6 =   4.45031289275240888144113950566e0;
const double Dopr853::b7 =   1.89151789931450038304281599044e0;
const double Dopr853::b8 =  -5.8012039600105847814672114227e0;
const double Dopr853::b9 =   3.1116436695781989440891606237e-1;
const double Dopr853::b10 = -1.52160949662516078556178806805e-1;
const double Dopr853::b11 =  2.01365400804030348374776537501e-1;
const double Dopr853::b12 =  4.47106157277725905176885569043e-2;
 
const double Dopr853::bhh1 = 0.244094488188976377952755905512e+00;
const double Dopr853::bhh2 = 0.733846688281611857341361741547e+00;
const double Dopr853::bhh3 = 0.220588235294117647058823529412e-01;

const double Dopr853::er1  =  0.1312004499419488073250102996e-01;
const double Dopr853::er6  = -0.1225156446376204440720569753e+01;
const double Dopr853::er7  = -0.4957589496572501915214079952e+00;
const double Dopr853::er8  =  0.1664377182454986536961530415e+01;
const double Dopr853::er9  = -0.3503288487499736816886487290e+00;
const double Dopr853::er10 =  0.3341791187130174790297318841e+00;
const double Dopr853::er11 =  0.8192320648511571246570742613e-01;
const double Dopr853::er12 = -0.2235530786388629525884427845e-01;
 
const double Dopr853::a21 =    5.26001519587677318785587544488e-2;
const double Dopr853::a31 =    1.97250569845378994544595329183e-2;
const double Dopr853::a32 =    5.91751709536136983633785987549e-2;
const double Dopr853::a41 =    2.95875854768068491816892993775e-2;
const double Dopr853::a43 =    8.87627564304205475450678981324e-2;
const double Dopr853::a51 =    2.41365134159266685502369798665e-1;
const double Dopr853::a53 =   -8.84549479328286085344864962717e-1;
const double Dopr853::a54 =    9.24834003261792003115737966543e-1;
const double Dopr853::a61 =    3.7037037037037037037037037037e-2;
const double Dopr853::a64 =    1.70828608729473871279604482173e-1;
const double Dopr853::a65 =    1.25467687566822425016691814123e-1;
const double Dopr853::a71 =    3.7109375e-2;
const double Dopr853::a74 =    1.70252211019544039314978060272e-1;
const double Dopr853::a75 =    6.02165389804559606850219397283e-2;
const double Dopr853::a76 =   -1.7578125e-2;
 
const double Dopr853::a81 =    3.70920001185047927108779319836e-2;
const double Dopr853::a84 =    1.70383925712239993810214054705e-1;
const double Dopr853::a85 =    1.07262030446373284651809199168e-1;
const double Dopr853::a86 =   -1.53194377486244017527936158236e-2;
const double Dopr853::a87 =    8.27378916381402288758473766002e-3;
const double Dopr853::a91 =    6.24110958716075717114429577812e-1;
const double Dopr853::a94 =   -3.36089262944694129406857109825e0;
const double Dopr853::a95 =   -8.68219346841726006818189891453e-1;
const double Dopr853::a96 =    2.75920996994467083049415600797e1;
const double Dopr853::a97 =    2.01540675504778934086186788979e1;
const double Dopr853::a98 =   -4.34898841810699588477366255144e1;
const double Dopr853::a101 =   4.77662536438264365890433908527e-1;
const double Dopr853::a104 =  -2.48811461997166764192642586468e0;
const double Dopr853::a105 =  -5.90290826836842996371446475743e-1;
const double Dopr853::a106 =   2.12300514481811942347288949897e1;
const double Dopr853::a107 =   1.52792336328824235832596922938e1;
const double Dopr853::a108 =  -3.32882109689848629194453265587e1;
const double Dopr853::a109 =  -2.03312017085086261358222928593e-2;
 
const double Dopr853::a111 =  -9.3714243008598732571704021658e-1;
const double Dopr853::a114 =   5.18637242884406370830023853209e0;
const double Dopr853::a115 =   1.09143734899672957818500254654e0;
const double Dopr853::a116 =  -8.14978701074692612513997267357e0;
const double Dopr853::a117 =  -1.85200656599969598641566180701e1;
const double Dopr853::a118 =   2.27394870993505042818970056734e1;
const double Dopr853::a119 =   2.49360555267965238987089396762e0;
const double Dopr853::a1110 = -3.0467644718982195003823669022e0;
const double Dopr853::a121 =   2.27331014751653820792359768449e0;
const double Dopr853::a124 =  -1.05344954667372501984066689879e1;
const double Dopr853::a125 =  -2.00087205822486249909675718444e0;
const double Dopr853::a126 =  -1.79589318631187989172765950534e1;
const double Dopr853::a127 =   2.79488845294199600508499808837e1;
const double Dopr853::a128 =  -2.85899827713502369474065508674e0;
const double Dopr853::a129 =  -8.87285693353062954433549289258e0;
const double Dopr853::a1210 =  1.23605671757943030647266201528e1;
const double Dopr853::a1211 =  6.43392746015763530355970484046e-1;
 
const double Dopr853::a141 =  5.61675022830479523392909219681e-2;
const double Dopr853::a147 =  2.53500210216624811088794765333e-1;
const double Dopr853::a148 = -2.46239037470802489917441475441e-1;
const double Dopr853::a149 = -1.24191423263816360469010140626e-1;
const double Dopr853::a1410 =  1.5329179827876569731206322685e-1;
const double Dopr853::a1411 =  8.20105229563468988491666602057e-3;
const double Dopr853::a1412 =  7.56789766054569976138603589584e-3;
const double Dopr853::a1413 = -8.298e-3;
 
const double Dopr853::a151 =  3.18346481635021405060768473261e-2;
const double Dopr853::a156 =  2.83009096723667755288322961402e-2;
const double Dopr853::a157 =  5.35419883074385676223797384372e-2;
const double Dopr853::a158 = -5.49237485713909884646569340306e-2;
const double Dopr853::a1511 = -1.08347328697249322858509316994e-4;
const double Dopr853::a1512 =  3.82571090835658412954920192323e-4;
const double Dopr853::a1513 = -3.40465008687404560802977114492e-4;
const double Dopr853::a1514 =  1.41312443674632500278074618366e-1;
const double Dopr853::a161 = -4.28896301583791923408573538692e-1;
const double Dopr853::a166 = -4.69762141536116384314449447206e0;
const double Dopr853::a167 =  7.68342119606259904184240953878e0;
const double Dopr853::a168 =  4.06898981839711007970213554331e0;
const double Dopr853::a169 =  3.56727187455281109270669543021e-1;
const double Dopr853::a1613 = -1.39902416515901462129418009734e-3;
const double Dopr853::a1614 =  2.9475147891527723389556272149e0;
const double Dopr853::a1615 = -9.15095847217987001081870187138e0;
 
const double Dopr853::d41  = -0.84289382761090128651353491142e+01;
const double Dopr853::d46  =  0.56671495351937776962531783590e+00;
const double Dopr853::d47  = -0.30689499459498916912797304727e+01;
const double Dopr853::d48  =  0.23846676565120698287728149680e+01;
const double Dopr853::d49  =  0.21170345824450282767155149946e+01;
const double Dopr853::d410 = -0.87139158377797299206789907490e+00;
const double Dopr853::d411 =  0.22404374302607882758541771650e+01;
const double Dopr853::d412 =  0.63157877876946881815570249290e+00;
const double Dopr853::d413 = -0.88990336451333310820698117400e-01;
const double Dopr853::d414 =  0.18148505520854727256656404962e+02;
const double Dopr853::d415 = -0.91946323924783554000451984436e+01;
const double Dopr853::d416 = -0.44360363875948939664310572000e+01;
 
const double Dopr853::d51  =  0.10427508642579134603413151009e+02;
const double Dopr853::d56  =  0.24228349177525818288430175319e+03;
const double Dopr853::d57  =  0.16520045171727028198505394887e+03;
const double Dopr853::d58  = -0.37454675472269020279518312152e+03;
const double Dopr853::d59  = -0.22113666853125306036270938578e+02;
const double Dopr853::d510 =  0.77334326684722638389603898808e+01;
const double Dopr853::d511 = -0.30674084731089398182061213626e+02;
const double Dopr853::d512 = -0.93321305264302278729567221706e+01;
const double Dopr853::d513 =  0.15697238121770843886131091075e+02;
const double Dopr853::d514 = -0.31139403219565177677282850411e+02;
const double Dopr853::d515 = -0.93529243588444783865713862664e+01;
const double Dopr853::d516 =  0.35816841486394083752465898540e+02;

const double Dopr853::d61 =  0.19985053242002433820987653617e+02;
const double Dopr853::d66 = -0.38703730874935176555105901742e+03;
const double Dopr853::d67 = -0.18917813819516756882830838328e+03;
const double Dopr853::d68 =  0.52780815920542364900561016686e+03;
const double Dopr853::d69 = -0.11573902539959630126141871134e+02;
const double Dopr853::d610 =  0.68812326946963000169666922661e+01;
const double Dopr853::d611 = -0.10006050966910838403183860980e+01;
const double Dopr853::d612 =  0.77771377980534432092869265740e+00;
const double Dopr853::d613 = -0.27782057523535084065932004339e+01;
const double Dopr853::d614 = -0.60196695231264120758267380846e+02;
const double Dopr853::d615 =  0.84320405506677161018159903784e+02;
const double Dopr853::d616 =  0.11992291136182789328035130030e+02;
 
const double Dopr853::d71  = -0.25693933462703749003312586129e+02;
const double Dopr853::d76  = -0.15418974869023643374053993627e+03;
const double Dopr853::d77  = -0.23152937917604549567536039109e+03;
const double Dopr853::d78  =  0.35763911791061412378285349910e+03;
const double Dopr853::d79  =  0.93405324183624310003907691704e+02;
const double Dopr853::d710 = -0.37458323136451633156875139351e+02;
const double Dopr853::d711 =  0.10409964950896230045147246184e+03;
const double Dopr853::d712 =  0.29840293426660503123344363579e+02;
const double Dopr853::d713 = -0.43533456590011143754432175058e+02;
const double Dopr853::d714 =  0.96324553959188282948394950600e+02;
const double Dopr853::d715 = -0.39177261675615439165231486172e+02;
const double Dopr853::d716 = -0.14972683625798562581422125276e+03;


void Dopr853::readParameterFile(std::ifstream &IN)
{
  IN >> startTime_;  // start time
  t_= startTime_;
  IN >> endTime_;    // end time
  
  IN >> printFlag_;  // output format
  IN >> numPrint_;   // number of time points printed to output
  
  IN >> h1_;   // max time step
  IN >> eps_;  // error tolerance
}

void Dopr853::simulate(void)
{ 
  // double tiny = 1e-30; // Using NR definition
  double tiny = 1e-9*eps_; // Caveat! Using new definition
  double  h, hNext, hDid;
  
  // Check that h1 and endTime - startTime are > 0
  //
  if (h1_ > 0.0 && (endTime_ - startTime_) > 0.0)
    h = h1_;
  else {//either h or (endTime-startTime) <=0
    std::cerr << "solverDopr853::simulate() - "
	      << "Wrong time borders or time step for simulation. "
	      << "No simulation performed.\n";
    exit(-1);
  }
  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  // Introduce the neighborhood
  //
  if (O_->numNeighborhood())
    O_->neighborhoodCreate(y_, t_);

  if (y_.size() && y_.size() != dydt_.size()) {
		std::cerr << "RK5Adaptive::simulate() "
							<< "resizing dydt_." << std::endl;
    dydt_.resize(y_.size(), y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
  // Create all vectors that will be needed here and by rkqs and rkck!
  //
  //Used here
  DataMatrix yScal( N() );
  //Used by rkqs
	DataMatrix yTemp( N() ),yErr( N() ),yErr2( N() );
  //used by rkck
  DataMatrix k2( N() ),k3( N() ),k4( N() ),
    k5( N() ),k6( N() ), k7( N() ), k8( N() ), k9( N() ), k10( N() ),
    yTempRkck( N() );
  //Resize each vector
  for (size_t i = 0; i < N(); i++) {
    yScal[i].resize(M());
    yErr[i].resize(M());
    yErr2[i].resize(M());
    yTemp[i].resize(M());
    k2[i].resize(M());
    k3[i].resize(M());
    k4[i].resize(M());
    k5[i].resize(M());
    k6[i].resize(M());
    k7[i].resize(M());
    k8[i].resize(M());
    k9[i].resize(M());
    k10[i].resize(M());
    yTempRkck[i].resize(M());
  }
  assert( y_.size()==yScal.size() && y_.size()==yErr.size() );
  assert( y_.size()==yErr2.size() && y_.size()==yTemp.size() );
  assert( y_.size()==k2.size() && y_.size()==k3.size() );
  assert( y_.size()==k4.size() && y_.size()==k5.size() );
  assert( y_.size()==k6.size() && y_.size()==k7.size() );
  assert( y_.size()==k8.size() && y_.size()==k9.size() );
  assert( y_.size()==k10.size() && y_.size()==yTempRkck.size() );
  
  // Initiate print times
  //
  double printTime = endTime_ + tiny;
  double printDeltaTime = endTime_ + 2.0 * tiny;
  //printFlag_=1;
  if (numPrint_ <= 0) //No printing
    printFlag_ = 0;
  else if (numPrint_ == 1) { // Print last point (default)
  }
  else if (numPrint_ == 2) { //Print first/last point
    printTime = startTime_ - tiny;
  } 
  else { //Print first/last points and spread the rest uniformly
    printTime = startTime_ - tiny;
    printDeltaTime = (endTime_ - startTime_) / ((double) (numPrint_ - 1));
  }

  // Go
  //////////////////////////////////////////////////////////////////////
  t_ = startTime_;
  numOk_ = numBad_ = 0;
  for (unsigned int nstp = 0;; nstp++) {
    if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
    // Update the derivatives
    O_->derivs(y_, dydt_);
    // Calculate 'scaling' for error measure
    for (size_t i = 0; i < N(); i++)
      for (size_t j = 0; j < M(); j++) {
        yScal[i][j] = fabs(y_[i][j]) + fabs(dydt_[i][j] * h) + tiny;
      }

    // Print if applicable 
    if (printFlag_ && t_ >= printTime) {
      printTime += printDeltaTime;
      print();
    }

    // Check if step is larger than max allowed
    // max step end is min of endTime_ and printTime
    double tMin = endTime_< printTime ? endTime_ : printTime;
    if (t_+h > tMin) h = tMin - t_;

    // Update
    rkqs(h, hDid, hNext, yScal, yTemp, yErr, yErr2, k2, k3, k4, k5, k6, 
	 k7, k8, k9, k10, yTempRkck);
    if (hDid == h) ++numOk_; else ++numBad_;
				
		// Make updates of rections that requires it
		O_->reactionUpdate(hDid, t_, y_);

    // Check for cell divisions, compartmental removal etc
        //DISABLED 042521
    //O_->compartmentChange(y_, dydt_, t_);
		
    // Rescale all temporary vectors as well
    if (y_.size() != yScal.size()) {
      yScal.resize(N(), yScal[0]);
      yTemp.resize(N(), yTemp[0]);
      yErr.resize(N(), yErr[0]);
      yErr2.resize(N(), yErr2[0]);
      k2.resize(N(), k2[0]);
      k3.resize(N(), k3[0]);
      k4.resize(N(), k4[0]);
      k5.resize(N(), k5[0]);
      k6.resize(N(), k6[0]);
      k7.resize(N(), k7[0]);
      k8.resize(N(), k8[0]);
      k9.resize(N(), k9[0]);
      k10.resize(N(), k10[0]);
      yTempRkck.resize(N(), yTempRkck[0]);
    }

    // Update neighborhood
    if (O_->numNeighborhood())
      O_->neighborhoodUpdate(y_, t_);
		
    // If the end t is passed return (print if applicable)
    if (t_ >= endTime_) {
      if (printFlag_) {
				// Update the derivatives
				O_->derivs(y_, dydt_);
				print();
      }
      std::cerr << "Simulation done.\n"; 
      return;
    }
    //Warn for small step sizes...
    //  if (fabs(hNext) <= hMin) {
    //        std::cerr << "Warning: Step size small (" << hNext
    //  		<< ") in rk5Adaptive::simulate at time " 
    //  		<< t << "\n"; /*exit(-1);*/ }
    h = hNext;
    //Do not take larger steps than h1
    if (h > h1_)
      h = h1_;
  }
}

#define SAFETY 0.85  //0.9
#define PGROW -0.1       //-0.2
#define PSHRNK -0.2    // -0.25
#define ERRCON 1.89e-4
void Dopr853::rkqs(double hTry, double &hDid, double &hNext,
											 std::vector< std::vector<double> > &yScal,
											 std::vector< std::vector<double> > &yTemp,
											 std::vector< std::vector<double> > &yErr,
											 std::vector< std::vector<double> > &yErr2,
											 std::vector< std::vector<double> > &k2,
											 std::vector< std::vector<double> > &k3,
											 std::vector< std::vector<double> > &k4,
											 std::vector< std::vector<double> > &k5,
											 std::vector< std::vector<double> > &k6,
											 std::vector< std::vector<double> > &k7,
											 std::vector< std::vector<double> > &k8,
											 std::vector< std::vector<double> > &k9,
											 std::vector< std::vector<double> > &k10,
											 std::vector< std::vector<double> > &yTempRkck)
{
  double errMax, h, hTemp, tNew, aux;
  h = hTry;
  for (;;) {   // infinite loop! ended by break statement
    rkck(h, yTemp, yErr, yErr2, k2, k3, k4, k5, k6, k7, k8, k9, k10, yTempRkck);

    errMax = 0.0;
    for (size_t i = 0; i < N(); i++)
      for (size_t j = 0; j < M(); j++) {
        aux = ((yErr2[i][j]*yErr2[i][j])/(sqrt(yErr[i][j]*yErr[i][j]+0.01*yErr2[i][j]*yErr2[i][j]))) / yScal[i][j];  // absolute value of the relative error
        if (aux > errMax)
          errMax = aux;  // save the largest error
      }
    errMax *= h;  // factor h not included when calculating the error, therefore: include it here
    errMax /= eps_;  // compare error to the tolerance level
    if (errMax <= 1.0) break;
    hTemp = SAFETY * h * pow(errMax, PSHRNK);
    if (h >= 0.0)
      h = hTemp > 0.1 * h ? hTemp : 0.1 * h;
    else
      h = hTemp > 0.1 * h ? 0.1 * h : hTemp;
    tNew = t_ + h;
    if (tNew == t_) { 
      std::cerr << "Warning stepsize underflow in solverDopr853::rkqs\n"; 
      exit(-1);
    }
  }
  if (errMax > ERRCON) hNext = SAFETY * h * pow(errMax, PGROW);
  else hNext = 5.0 * h;
  t_ += (hDid = h);
  
  for (size_t i = 0; i < N(); i++)
    for (size_t j = 0; j < M(); j++)
      y_[i][j] = yTemp[i][j];
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

void Dopr853::rkck(double h,
		       std::vector< std::vector<double> > &yOut,
		       std::vector< std::vector<double> > &yErr,
		       std::vector< std::vector<double> > &yErr2,
		       std::vector< std::vector<double> > &k2,
		       std::vector< std::vector<double> > &k3,
		       std::vector< std::vector<double> > &k4,
		       std::vector< std::vector<double> > &k5,
		       std::vector< std::vector<double> > &k6,
		       std::vector< std::vector<double> > &k7,
		       std::vector< std::vector<double> > &k8,
		       std::vector< std::vector<double> > &k9,
		       std::vector< std::vector<double> > &k10,
		       std::vector< std::vector<double> > &yTempRkck ) {
  
  double n = N();
  double m = M();
	
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+a21*h*dydt_[i][j];
  
  O_->derivs(yTempRkck, k2); // t + a2 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*(a31*dydt_[i][j]+a32*k2[i][j]);
  
  O_->derivs(yTempRkck, k3); // t + a3 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*(a41*dydt_[i][j]+a43*k3[i][j]);
  
  O_->derivs(yTempRkck, k4); // t + a4 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*
	(a51*dydt_[i][j]+a53*k3[i][j]+a54*k4[i][j]);
  
  O_->derivs(yTempRkck, k5); // t + a5 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*
	(a61*dydt_[i][j]+a64*k4[i][j]+a65*k5[i][j]);
  
  O_->derivs(yTempRkck, k6); // t + a6 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*
	(a71*dydt_[i][j]+a74*k4[i][j]+a75*k5[i][j]+a76*k6[i][j]);

  O_->derivs(yTempRkck, k7); // t + a7 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*
	(a81*dydt_[i][j]+a84*k4[i][j]+a85*k5[i][j]+a86*k6[i][j]+
	 a87*k7[i][j]);

  O_->derivs(yTempRkck, k8); // t + a8 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*
	(a91*dydt_[i][j]+a94*k4[i][j]+a95*k5[i][j]+a96*k6[i][j]+
	 a97*k7[i][j]+a98*k8[i][j]);

  O_->derivs(yTempRkck, k9); // t + a9 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*
	(a101*dydt_[i][j]+a104*k4[i][j]+a105*k5[i][j]+a106*k6[i][j]+
	 a107*k7[i][j]+a108*k8[i][j]+a109*k9[i][j]);

  O_->derivs(yTempRkck, k10); // t + a10 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*
	(a111*dydt_[i][j]+a114*k4[i][j]+a115*k5[i][j]+a116*k6[i][j]
	 +a117*k7[i][j]+a118*k8[i][j]+a119*k9[i][j]+a1110*k10[i][j]);
  
  O_->derivs(yTempRkck, k2); // t + a11 * h
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+h*
	(a121*dydt_[i][j]+a124*k4[i][j]+a125*k5[i][j]+a126*k6[i][j]
	 +a127*k7[i][j]+a128*k8[i][j]+a129*k9[i][j]+a1210*k10[i][j]
	 +a1211*k2[i][j]);
  
  O_->derivs(yTempRkck, k3); // t + h
  for (size_t i = 0; i < n; ++i){
    for (size_t j = 0; j < m; ++j){
      k4[i][j]=b1*dydt_[i][j]+b6*k6[i][j]+b7*k7[i][j]+b8*k8[i][j]+b9*k9[i][j]+b10*k10[i][j]+b11*k2[i][j]+b12*k3[i][j];
      yOut[i][j]=y_[i][j]+h*k4[i][j];
    }
  }

  for (size_t i = 0; i < n; ++i){
    for (size_t j = 0; j < m; ++j){
      yErr[i][j]=k4[i][j]-bhh1*dydt_[i][j]-bhh2*k9[i][j]-bhh3*k3[i][j];   // the 3rd order error
      yErr2[i][j]=er1*dydt_[i][j]+er6*k6[i][j]+er7*k7[i][j]+er8*k8[i][j]+er9*k9[i][j]+
	er10*k10[i][j]+er11*k2[i][j]+er12*k3[i][j];  // the 5th order error
    }
  }
}

///////////////////////////////////////////////////////////////////
//                                                               //
// RK5Adaptive: 5th order Runge-kutta with adaptive step-length  //
//                                                               //
///////////////////////////////////////////////////////////////////


RK5Adaptive::RK5Adaptive(Organism *O,std::ifstream &IN)
  :BaseSolver(O,IN)
{
	readParameterFile(IN);
}

void RK5Adaptive::readParameterFile(std::ifstream &IN)
{
  IN >> startTime_;
  t_= startTime_;
  IN >> endTime_;
  
  IN >> printFlag_;
  IN >> numPrint_;
  
  IN >> h1_;
  IN >> eps_;
}

void RK5Adaptive::simulate(void)
{ 
  // double tiny = 1e-30; // Using NR definition
  double tiny = 1e-9*eps_; // Caveat! Using new definition
  double  h, hNext, hDid;
  
  // Check that h1 and endTime - startTime are > 0
  //
  if (h1_ > 0.0 && (endTime_ - startTime_) > 0.0)
    h = h1_;
  else {//either h or (endTime-startTime) <=0
    std::cerr << "solverRk5Adaptive::simulate() - "
	      << "Wrong time borders or time step for simulation. "
	      << "No simulation performed.\n";
    exit(-1);
  }
  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  // Introduce the neighborhood
  //
  if (O_->numNeighborhood())
    O_->neighborhoodCreate(y_, t_);
  if (y_.size() && y_.size() != dydt_.size()) {
    dydt_.resize(y_.size(), y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
  // Create all vectors that will be needed here and by rkqs and rkck!
  //
  //Used here
  std::vector< std::vector<double> > yScal( N() );
  //Used by rkqs
  std::vector< std::vector<double> > yTemp( N() ),yErr( N() );
  //used by rkck
  std::vector< std::vector<double> > ak2( N() ),ak3( N() ),ak4( N() ),
    ak5( N() ),ak6( N() ),yTempRkck( N() );
  //Resize each vector
  for (size_t i = 0; i < N(); i++) {
    yScal[i].resize(M());
    yErr[i].resize(M());
    yTemp[i].resize(M());
    ak2[i].resize(M());
    ak3[i].resize(M());
    ak4[i].resize(M());
    ak5[i].resize(M());
    ak6[i].resize(M());
    yTempRkck[i].resize(M());
  }
  assert( y_.size()==yScal.size() && y_.size()==yErr.size() );
  assert( y_.size()==yTemp.size() && y_.size()==ak2.size() );
  assert( y_.size()==ak3.size() && y_.size()==ak4.size() );
  assert( y_.size()==ak5.size() && y_.size()==ak6.size() );
  assert( y_.size()==yTempRkck.size() );
  
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
  double printTime = endTime_ + tiny;
  double printDeltaTime = endTime_ + 2.0 * tiny;
  if (numPrint_ <= 0) //No printing
    printFlag_ = 0;
  else if (numPrint_ == 1) { // Print last point (default)
  }
  else if (numPrint_ == 2) { //Print first/last point
    printTime = startTime_ - tiny;
  } 
  else { //Print first/last points and spread the rest uniformly
    printTime = startTime_ - tiny;
    printDeltaTime = (endTime_ - startTime_) / ((double) (numPrint_ - 1));
  }

  // Go
  //////////////////////////////////////////////////////////////////////
  t_ = startTime_;
  numOk_ = numBad_ = 0;
  for (unsigned int nstp = 0;; nstp++) {
    if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
    // Update the derivatives
    O_->derivs(y_, dydt_);
    // Calculate 'scaling' for error measure
    for (size_t i = 0; i < N(); i++)
      for (size_t j = 0; j < M(); j++) {
        yScal[i][j] = fabs(y_[i][j]) + fabs(dydt_[i][j] * h) + tiny;
      }

    // Print if applicable 
    if (printFlag_ && t_ >= printTime) {
      printTime += printDeltaTime;
      print();
    }

    // Check if step is larger than max allowed
    // max step end is min of endTime_ and printTime
    double tMin = endTime_< printTime ? endTime_ : printTime;
    if (t_+h > tMin) h = tMin - t_;

    // Update
    rkqs(h, hDid, hNext, yScal, yTemp, yErr, ak2, ak3, ak4, ak5, ak6, 
				 yTempRkck);
    if (hDid == h) ++numOk_; else ++numBad_;
				
		// Make updates of rections that requires it
		O_->reactionUpdate(hDid, t_, y_);

    // Check for cell divisions, compartmental removal etc
        //DISABLED 042521
    //O_->compartmentChange(y_, dydt_, t_);
		
    // Rescale all temporary vectors as well
    if (y_.size() != yScal.size()) {
      yScal.resize(N(), yScal[0]);
      yTemp.resize(N(), yTemp[0]);
      yErr.resize(N(), yErr[0]);
      ak2.resize(N(), ak2[0]);
      ak3.resize(N(), ak3[0]);
      ak4.resize(N(), ak4[0]);
      ak5.resize(N(), ak5[0]);
      ak6.resize(N(), ak6[0]);
      yTempRkck.resize(N(), yTempRkck[0]);
    }
		
    // Update neighborhood
    if (O_->numNeighborhood())
      O_->neighborhoodUpdate(y_, t_);
		
    // If the end t is passed return (print if applicable)
    if (t_ >= endTime_) {
      if (printFlag_) {
				// Update the derivatives
				O_->derivs(y_, dydt_);
				print();
      }
      std::cerr << "Simulation done.\n"; 
      return;
    }
    //Warn for small step sizes...
    //  if (fabs(hNext) <= hMin) {
    //        std::cerr << "Warning: Step size small (" << hNext
    //  		<< ") in rk5Adaptive::simulate at time " 
    //  		<< t << "\n"; /*exit(-1);*/ }
    h = hNext;
    //Do not take larger steps than h1
    if (h > h1_)
      h = h1_;
  }
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void RK5Adaptive::rkqs(double hTry, double &hDid, double &hNext,
		std::vector< std::vector<double> > &yScal,
		std::vector< std::vector<double> > &yTemp,
		std::vector< std::vector<double> > &yErr,
		std::vector< std::vector<double> > &ak2,
		std::vector< std::vector<double> > &ak3,
		std::vector< std::vector<double> > &ak4,
		std::vector< std::vector<double> > &ak5,
		std::vector< std::vector<double> > &ak6,
		std::vector< std::vector<double> > &yTempRkck)
{
	double errMax, h, hTemp, tNew, aux;
	h = hTry;
	double nanProtection = 10000.*eps_; // 2.*eps_; // must be larger than 1/eps_
	for (;;) {
		rkck(h, yTemp, yErr, ak2, ak3, ak4, ak5, ak6, yTempRkck);
		errMax = 0.0;
		for (size_t i = 0; i < N(); i++)
			for (size_t j = 0; j < M(); j++) {
				aux = fabs(yErr[i][j] / yScal[i][j]);
				if (aux > errMax) // Warning: if aux is nan then this will evaluate as false.
					errMax = aux;
				else if (std::isnan(aux)){ // if aux==nan, ensure that we reduce step size.
					std::cerr << "NAN forced stepsize decrease in rungeKutta.cc RK5Adaptive::rkqs\n";
					if (errMax <  nanProtection)
						errMax = nanProtection;
				}
			}
		errMax /= eps_;
		if (errMax <= 1.0) break;
		hTemp = SAFETY * h * pow(errMax, PSHRNK);
		if (h >= 0.0)
			h = hTemp > 0.1 * h ? hTemp : 0.1 * h;
		else
			h = hTemp > 0.1 * h ? 0.1 * h : hTemp;
		tNew = t_ + h;
		if (tNew == t_) { 
			std::cerr << "Warning stepsize underflow in solverRk5Adaptive::rkqs\n"; 
			exit(-1);
		}
	}
	if (errMax > ERRCON) hNext = SAFETY * h * pow(errMax, PGROW);
	else hNext = 5.0 * h;
	t_ += (hDid = h);

	for (size_t i = 0; i < N(); i++)
		for (size_t j = 0; j < M(); j++)
			y_[i][j] = yTemp[i][j];
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

void RK5Adaptive::
updateMechanicsUntilConvergence(double cEps,double cTmax,
				std::vector< std::vector<double> > &yScal,
				std::vector< std::vector<double> > &yTemp,
				std::vector< std::vector<double> > &yErr,
				std::vector< std::vector<double> > &ak2,
				std::vector< std::vector<double> > &ak3,
				std::vector< std::vector<double> > &ak4,
				std::vector< std::vector<double> > &ak5,
				std::vector< std::vector<double> > &ak6,
				std::vector< std::vector<double> > &yTempRkck,
				std::vector<size_t> &mechEq,
				std::vector<size_t> &mechPos ) 
{  
  //double tiny=1e-30;//Using NR definition
  double tiny=1e-9*eps_;//Caveat! Using new definition
  double h=h1_,hNext,hDid;

  // Go
  //////////////////////////////////////////////////////////////////////
  double t=0.0;
  //numOk_ = numBad_ = 0;
  for( unsigned int nstp=0 ; ; nstp++ ) {
		if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
    //Update the derivatives
    O_->derivsMechanical(y_,dydt_,mechEq,mechPos);
    
    //Calculate 'scaling' for error measure
    for( size_t i=0 ; i<N() ; i++ )
      for( size_t k=0 ; k<mechPos.size() ; k++ ) {
				size_t j=mechPos[k];
        yScal[i][j]=fabs(y_[i][j])+fabs(dydt_[i][j]*h)+tiny;
      }
    
    //Update
    rkqsMechanical(t,h,hDid,hNext,yScal,
		   yTemp,yErr,ak2,ak3,ak4,ak5,ak6,yTempRkck,
		   mechEq,mechPos);
    //if (hDid == h) ++numOk_; else ++numBad_;
    
		// Make updates of rections that requires it
		O_->reactionUpdate(hDid, t_, y_);
    
    //Update neighborhood
    if( O_->numNeighborhood() )
      O_->neighborhoodUpdate(y_,t_);

    //If the positional variables have converged 
    //or maximal T is reached, return
    double maxDydt=0.0;
    for( size_t i=0 ; i<N() ; i++ )
      for( size_t k=0 ; k<mechPos.size() ; k++ ) {
	size_t j=mechPos[k];
	if( std::fabs(dydt_[i][j])>maxDydt )
	  maxDydt = std::fabs(dydt_[i][j]); 
      }
    if( maxDydt<cEps ) {
      //std::cerr << "Rk5Adaptive::updateMechanicsUntilConvergance "
      //	<< "t=" << t << " maxDydt=" << maxDydt << "\n";
      return;
    }
    else if ( t>=cTmax ) {
      std::cerr << "Rk5Adaptive::updateMechanicsUntilConvergance\n"
		<< "Warning, final time " << cTmax << " reached with"
		<< " Max dy/dt " << maxDydt << ".\n"; 
      return;
    }
    //Warn for small step sizes...
    //  if (fabs(hNext) <= hMin) {
    //        std::cerr << "Warning: Step size small (" << hNext
    //  		<< ") in rk5Adaptive::simulate at time " 
    //  		<< t << "\n"; /*exit(-1);*/ }
    h=hNext;
    //Do not take larger steps than h1
    if( h>h1_ )
      h=h1_;
  }  
}


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
void RK5Adaptive::
rkqsMechanical(double &t,double hTry,double &hDid,double &hNext,
	       std::vector< std::vector<double> > &yScal,
	       std::vector< std::vector<double> > &yTemp,
	       std::vector< std::vector<double> > &yErr,
	       std::vector< std::vector<double> > &ak2,
	       std::vector< std::vector<double> > &ak3,
	       std::vector< std::vector<double> > &ak4,
	       std::vector< std::vector<double> > &ak5,
	       std::vector< std::vector<double> > &ak6,
	       std::vector< std::vector<double> > &yTempRkck,
	       std::vector<size_t> &mechEq,std::vector<size_t> &posVar) 
{  
  double errMax,h,hTemp,tNew,aux;
  h=hTry;
  for (;;) {
    rkckMechanical(h,yTemp,yErr,ak2,ak3,ak4,ak5,ak6,yTempRkck,
		   mechEq,posVar);
    errMax=0.;
    for(size_t i=0 ; i<N() ; i++ )
      for( size_t k=0 ; k<posVar.size() ; k++ ) {
	size_t j=posVar[k];
        aux = fabs(yErr[i][j]/yScal[i][j]);
        if( aux>errMax )
          errMax = aux;
      }
    errMax /= eps_;
    if (errMax <= 1.0 ) break;
    hTemp=SAFETY*h*pow(errMax,PSHRNK);
    if( h >= 0. )
      h = hTemp > 0.1*h ? hTemp : 0.1*h;
    else
      h = hTemp > 0.1*h ? 0.1*h : hTemp;
    tNew=t+h;
    if( tNew==t ) { 
      std::cerr << "Warning stepsize underflow in Rk5Adaptive::rkqs\n"; 
      exit(-1); }
  }
  if (errMax > ERRCON) hNext=SAFETY*h*pow(errMax,PGROW);
  else hNext=5.0*h;
  t += (hDid=h);
  
  for(size_t i=0 ; i<N() ; i++ )
    for( size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      y_[i][j]=yTemp[i][j];
    }
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

void RK5Adaptive::rkck(double h,
		       std::vector< std::vector<double> > &yOut,
		       std::vector< std::vector<double> > &yErr,
		       std::vector< std::vector<double> > &ak2,
		       std::vector< std::vector<double> > &ak3,
		       std::vector< std::vector<double> > &ak4,
		       std::vector< std::vector<double> > &ak5,
		       std::vector< std::vector<double> > &ak6,
		       std::vector< std::vector<double> > &yTempRkck ) {
  
  //static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875;
  static double b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  double n = N();
	double m = M();
	
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+b21*h*dydt_[i][j];
  
  O_->derivs(yTempRkck, ak2); // t + a2h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yTempRkck[i][j]=y_[i][j]+h*(b31*dydt_[i][j]+b32*ak2[i][j]);
  
  O_->derivs(yTempRkck, ak3); // t + a3h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yTempRkck[i][j]=y_[i][j]+h*(b41*dydt_[i][j]+b42*ak2[i][j]+b43*ak3[i][j]);
  
  O_->derivs(yTempRkck, ak4); // t + a4 * h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yTempRkck[i][j]=y_[i][j]+h*
				(b51*dydt_[i][j]+b52*ak2[i][j]+b53*ak3[i][j]+b54*ak4[i][j]);
  
  O_->derivs(yTempRkck, ak5); // t + a5 * h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yTempRkck[i][j]=y_[i][j]+h*
				(b61*dydt_[i][j]+b62*ak2[i][j]+b63*ak3[i][j]+b64*ak4[i][j]+
				 b65*ak5[i][j]);
  
  O_->derivs(yTempRkck, ak6); // t + a6 * h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j) {
			yOut[i][j]=y_[i][j]+h*
				(c1*dydt_[i][j]+c3*ak3[i][j]+c4*ak4[i][j]+c6*ak6[i][j]);
		}
  
	for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j) 
      yErr[i][j]=h*
        (dc1*dydt_[i][j]+dc3*ak3[i][j]+dc4*ak4[i][j]+
         dc5*ak5[i][j]+dc6*ak6[i][j]);

}

void RK5Adaptive::
rkckMechanical(double h,
	       std::vector< std::vector<double> > &yOut,
	       std::vector< std::vector<double> > &yErr,
	       std::vector< std::vector<double> > &ak2,
	       std::vector< std::vector<double> > &ak3,
	       std::vector< std::vector<double> > &ak4,
	       std::vector< std::vector<double> > &ak5,
	       std::vector< std::vector<double> > &ak6,
	       std::vector< std::vector<double> > &yTempRkck,
	       std::vector<size_t> &mechEq,std::vector<size_t> &posVar ) 
{  
  //static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875;
  static double b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+b21*h*dydt_[i][j];
    }
  O_->derivsMechanical(yTempRkck,ak2,mechEq,posVar);//t+a2h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+h*(b31*dydt_[i][j]+b32*ak2[i][j]);
    }
  O_->derivsMechanical(yTempRkck,ak3,mechEq,posVar);//t+a3h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+h*(b41*dydt_[i][j]+b42*ak2[i][j]+b43*ak3[i][j]);
    }
  O_->derivsMechanical(yTempRkck,ak4,mechEq,posVar);//t+a4*h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+h*
        (b51*dydt_[i][j]+b52*ak2[i][j]+b53*ak3[i][j]+b54*ak4[i][j]);
    }
  O_->derivsMechanical(yTempRkck,ak5,mechEq,posVar);//t+a5*h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+h*
        (b61*dydt_[i][j]+b62*ak2[i][j]+b63*ak3[i][j]+b64*ak4[i][j]+
         b65*ak5[i][j]);
    }  
  O_->derivsMechanical(yTempRkck,ak6,mechEq,posVar);//t+a6*h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yOut[i][j]=y_[i][j]+h*
        (c1*dydt_[i][j]+c3*ak3[i][j]+c4*ak4[i][j]+c6*ak6[i][j]);
    }
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yErr[i][j]=h*
        (dc1*dydt_[i][j]+dc3*ak3[i][j]+dc4*ak4[i][j]+
         dc5*ak5[i][j]+dc6*ak6[i][j]);
    }
}

///////////////////////////////////////////////////////////////////
//                                                               //
// RK5AdaptiveEquilibrium: 5th order Runge-kutta with adaptive   //
// step-length, simulates to an equilibrium condition            //
//                                                               //
///////////////////////////////////////////////////////////////////


RK5AdaptiveEquilibrium::RK5AdaptiveEquilibrium(Organism *O,std::ifstream &IN)
  :BaseSolver(O,IN)
{
	readParameterFile(IN);
}

void RK5AdaptiveEquilibrium::readParameterFile(std::ifstream &IN)
{
  IN >> startTime_;
  t_= startTime_;
  IN >> endTime_;
  
  IN >> printFlag_;
  IN >> printInterval_;
  numPrint_ = 0;
  
  IN >> h1_;
  IN >> eps_;
  IN >> threshold_;
}

void RK5AdaptiveEquilibrium::simulate(void)
{ 
  // double tiny = 1e-30; // Using NR definition
  double tiny = 1e-9*eps_; // Caveat! Using new definition
  double  h, hNext, hDid;
  
  // Check that h1 and endTime - startTime are > 0
  //
  if (h1_ > 0.0 && (endTime_ - startTime_) > 0.0)
    h = h1_;
  else {//either h or (endTime-startTime) <=0
    std::cerr << "solverRk5Adaptive::simulate() - "
	      << "Wrong time borders or time step for simulation. "
	      << "No simulation performed.\n";
    exit(-1);
  }
  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  // Introduce the neighborhood
  //
  if (O_->numNeighborhood())
    O_->neighborhoodCreate(y_, t_);
  if (y_.size() && y_.size() != dydt_.size()) {
    dydt_.resize(y_.size(), y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
  // Create all vectors that will be needed here and by rkqs and rkck!
  //
  //Used here
  std::vector< std::vector<double> > yScal( N() );
  //Used by rkqs
  std::vector< std::vector<double> > yTemp( N() ),yErr( N() );
  //used by rkck
  std::vector< std::vector<double> > ak2( N() ),ak3( N() ),ak4( N() ),
    ak5( N() ),ak6( N() ),yTempRkck( N() );
  //Resize each vector
  for (size_t i = 0; i < N(); i++) {
    yScal[i].resize(M());
    yErr[i].resize(M());
    yTemp[i].resize(M());
    ak2[i].resize(M());
    ak3[i].resize(M());
    ak4[i].resize(M());
    ak5[i].resize(M());
    ak6[i].resize(M());
    yTempRkck[i].resize(M());
  }
  assert( y_.size()==yScal.size() && y_.size()==yErr.size() );
  assert( y_.size()==yTemp.size() && y_.size()==ak2.size() );
  assert( y_.size()==ak3.size() && y_.size()==ak4.size() );
  assert( y_.size()==ak5.size() && y_.size()==ak6.size() );
  assert( y_.size()==yTempRkck.size() );
  
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
  double printTime;
  double printDeltaTime=printInterval_;
  if( printInterval_ == 0 )// Print only last point
    printTime=endTime_ + tiny;
  else {//Print first/last points and spread the rest uniformly
    printTime=startTime_-tiny;

  }

  // Go
  //////////////////////////////////////////////////////////////////////
  double maxDeriv = 1000.0;
  double tempDeriv = 1000.0;
  t_ = startTime_;
  numOk_ = numBad_ = 0;
  for (unsigned int nstp = 0;; nstp++) {
    if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
    // Update the derivatives
    O_->derivs(y_, dydt_);
    // Calculate 'scaling' for error measure
    for (size_t i = 0; i < N(); i++)
      for (size_t j = 0; j < M(); j++) {
        yScal[i][j] = fabs(y_[i][j]) + fabs(dydt_[i][j] * h) + tiny;
      }

    //Print if applicable 
    if( printFlag_ && t_ >= printTime ) {
      printTime += printDeltaTime;
      print();
      numPrint_ +=1;
    }

    // Check if step is larger than max allowed
    // max step end is min of endTime_ and printTime
    double tMin = endTime_< printTime ? endTime_ : printTime;
    if (t_+h > tMin) h = tMin - t_;

    // Update
    rkqs(h, hDid, hNext, yScal, yTemp, yErr, ak2, ak3, ak4, ak5, ak6, 
				 yTempRkck);
    if (hDid == h) ++numOk_; else ++numBad_;
				
		// Make updates of rections that requires it
		O_->reactionUpdate(hDid, t_, y_);

    // Check for cell divisions, compartmental removal etc
        //DISABLED 042521
    //O_->compartmentChange(y_, dydt_, t_);
		
    // Rescale all temporary vectors as well
    if (y_.size() != yScal.size()) {
      yScal.resize(N(), yScal[0]);
      yTemp.resize(N(), yTemp[0]);
      yErr.resize(N(), yErr[0]);
      ak2.resize(N(), ak2[0]);
      ak3.resize(N(), ak3[0]);
      ak4.resize(N(), ak4[0]);
      ak5.resize(N(), ak5[0]);
      ak6.resize(N(), ak6[0]);
      yTempRkck.resize(N(), yTempRkck[0]);
    }
		
    // Update neighborhood
    if (O_->numNeighborhood())
      O_->neighborhoodUpdate(y_, t_);

    // find the maxDeriv
    maxDeriv = 0;
    for (size_t i = 0; i < N(); i++){
      for (size_t j = 0; j < M(); j++) {
	tempDeriv = fabs(dydt_[i][j]);
	if(tempDeriv > maxDeriv){
	  maxDeriv = tempDeriv;	}
      }
    }
		
    // If the end t is passed, or maxDeriv is lower
    // than the threshold, return (print if applicable)
    if (t_ >= endTime_ || maxDeriv < threshold_) {
      if (printFlag_) {
				// Update the derivatives
				O_->derivs(y_, dydt_);
				print();
      }
      std::cerr << "Simulation done.\n"; 
      return;
    }
    //Warn for small step sizes...
    //  if (fabs(hNext) <= hMin) {
    //        std::cerr << "Warning: Step size small (" << hNext
    //  		<< ") in rk5Adaptive::simulate at time " 
    //  		<< t << "\n"; /*exit(-1);*/ }
    h = hNext;
    //Do not take larger steps than h1
    if (h > h1_)
      h = h1_;
  }
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void RK5AdaptiveEquilibrium::rkqs(double hTry, double &hDid, double &hNext,
											 std::vector< std::vector<double> > &yScal,
											 std::vector< std::vector<double> > &yTemp,
											 std::vector< std::vector<double> > &yErr,
											 std::vector< std::vector<double> > &ak2,
											 std::vector< std::vector<double> > &ak3,
											 std::vector< std::vector<double> > &ak4,
											 std::vector< std::vector<double> > &ak5,
											 std::vector< std::vector<double> > &ak6,
											 std::vector< std::vector<double> > &yTempRkck)
{
  double errMax, h, hTemp, tNew, aux;
  h = hTry;
  for (;;) {
    rkck(h, yTemp, yErr, ak2, ak3, ak4, ak5, ak6, yTempRkck);
    errMax = 0.0;
    for (size_t i = 0; i < N(); i++)
      for (size_t j = 0; j < M(); j++) {
        aux = fabs(yErr[i][j] / yScal[i][j]);
        if (aux > errMax)
          errMax = aux;
      }
    errMax /= eps_;
    if (errMax <= 1.0) break;
    hTemp = SAFETY * h * pow(errMax, PSHRNK);
    if (h >= 0.0)
      h = hTemp > 0.1 * h ? hTemp : 0.1 * h;
    else
      h = hTemp > 0.1 * h ? 0.1 * h : hTemp;
    tNew = t_ + h;
    if (tNew == t_) { 
      std::cerr << "Warning stepsize underflow in solverRk5Adaptive::rkqs\n"; 
      exit(-1);
    }
  }
  if (errMax > ERRCON) hNext = SAFETY * h * pow(errMax, PGROW);
  else hNext = 5.0 * h;
  t_ += (hDid = h);
  
  for (size_t i = 0; i < N(); i++)
    for (size_t j = 0; j < M(); j++)
      y_[i][j] = yTemp[i][j];
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

void RK5AdaptiveEquilibrium::
updateMechanicsUntilConvergence(double cEps,double cTmax,
				std::vector< std::vector<double> > &yScal,
				std::vector< std::vector<double> > &yTemp,
				std::vector< std::vector<double> > &yErr,
				std::vector< std::vector<double> > &ak2,
				std::vector< std::vector<double> > &ak3,
				std::vector< std::vector<double> > &ak4,
				std::vector< std::vector<double> > &ak5,
				std::vector< std::vector<double> > &ak6,
				std::vector< std::vector<double> > &yTempRkck,
				std::vector<size_t> &mechEq,
				std::vector<size_t> &mechPos ) 
{  
  //double tiny=1e-30;//Using NR definition
  double tiny=1e-9*eps_;//Caveat! Using new definition
  double h=h1_,hNext,hDid;

  // Go
  //////////////////////////////////////////////////////////////////////
  double t=0.0;
  //numOk_ = numBad_ = 0;
  for( unsigned int nstp=0 ; ; nstp++ ) {
		if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
    //Update the derivatives
    O_->derivsMechanical(y_,dydt_,mechEq,mechPos);
    
    //Calculate 'scaling' for error measure
    for( size_t i=0 ; i<N() ; i++ )
      for( size_t k=0 ; k<mechPos.size() ; k++ ) {
				size_t j=mechPos[k];
        yScal[i][j]=fabs(y_[i][j])+fabs(dydt_[i][j]*h)+tiny;
      }
    
    //Update
    rkqsMechanical(t,h,hDid,hNext,yScal,
		   yTemp,yErr,ak2,ak3,ak4,ak5,ak6,yTempRkck,
		   mechEq,mechPos);
    //if (hDid == h) ++numOk_; else ++numBad_;
    
		// Make updates of rections that requires it
		O_->reactionUpdate(hDid, t_, y_);
    
    //Update neighborhood
    if( O_->numNeighborhood() )
      O_->neighborhoodUpdate(y_,t_);

    //If the positional variables have converged 
    //or maximal T is reached, return
    double maxDydt=0.0;
    for( size_t i=0 ; i<N() ; i++ )
      for( size_t k=0 ; k<mechPos.size() ; k++ ) {
	size_t j=mechPos[k];
	if( std::fabs(dydt_[i][j])>maxDydt )
	  maxDydt = std::fabs(dydt_[i][j]); 
      }
    if( maxDydt<cEps ) {
      //std::cerr << "Rk5Adaptive::updateMechanicsUntilConvergance "
      //	<< "t=" << t << " maxDydt=" << maxDydt << "\n";
      return;
    }
    else if ( t>=cTmax ) {
      std::cerr << "Rk5Adaptive::updateMechanicsUntilConvergance\n"
		<< "Warning, final time " << cTmax << " reached with"
		<< " Max dy/dt " << maxDydt << ".\n"; 
      return;
    }
    //Warn for small step sizes...
    //  if (fabs(hNext) <= hMin) {
    //        std::cerr << "Warning: Step size small (" << hNext
    //  		<< ") in rk5Adaptive::simulate at time " 
    //  		<< t << "\n"; /*exit(-1);*/ }
    h=hNext;
    //Do not take larger steps than h1
    if( h>h1_ )
      h=h1_;
  }  
}


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
void RK5AdaptiveEquilibrium::
rkqsMechanical(double &t,double hTry,double &hDid,double &hNext,
	       std::vector< std::vector<double> > &yScal,
	       std::vector< std::vector<double> > &yTemp,
	       std::vector< std::vector<double> > &yErr,
	       std::vector< std::vector<double> > &ak2,
	       std::vector< std::vector<double> > &ak3,
	       std::vector< std::vector<double> > &ak4,
	       std::vector< std::vector<double> > &ak5,
	       std::vector< std::vector<double> > &ak6,
	       std::vector< std::vector<double> > &yTempRkck,
	       std::vector<size_t> &mechEq,std::vector<size_t> &posVar) 
{  
  double errMax,h,hTemp,tNew,aux;
  h=hTry;
  for (;;) {
    rkckMechanical(h,yTemp,yErr,ak2,ak3,ak4,ak5,ak6,yTempRkck,
		   mechEq,posVar);
    errMax=0.;
    for(size_t i=0 ; i<N() ; i++ )
      for( size_t k=0 ; k<posVar.size() ; k++ ) {
	size_t j=posVar[k];
        aux = fabs(yErr[i][j]/yScal[i][j]);
        if( aux>errMax )
          errMax = aux;
      }
    errMax /= eps_;
    if (errMax <= 1.0 ) break;
    hTemp=SAFETY*h*pow(errMax,PSHRNK);
    if( h >= 0. )
      h = hTemp > 0.1*h ? hTemp : 0.1*h;
    else
      h = hTemp > 0.1*h ? 0.1*h : hTemp;
    tNew=t+h;
    if( tNew==t ) { 
      std::cerr << "Warning stepsize underflow in Rk5Adaptive::rkqs\n"; 
      exit(-1); }
  }
  if (errMax > ERRCON) hNext=SAFETY*h*pow(errMax,PGROW);
  else hNext=5.0*h;
  t += (hDid=h);
  
  for(size_t i=0 ; i<N() ; i++ )
    for( size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      y_[i][j]=yTemp[i][j];
    }
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

void RK5AdaptiveEquilibrium::rkck(double h,
		       std::vector< std::vector<double> > &yOut,
		       std::vector< std::vector<double> > &yErr,
		       std::vector< std::vector<double> > &ak2,
		       std::vector< std::vector<double> > &ak3,
		       std::vector< std::vector<double> > &ak4,
		       std::vector< std::vector<double> > &ak5,
		       std::vector< std::vector<double> > &ak6,
		       std::vector< std::vector<double> > &yTempRkck ) {
  
  //static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875;
  static double b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  double n = N();
	double m = M();
	
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yTempRkck[i][j]=y_[i][j]+b21*h*dydt_[i][j];
  
  O_->derivs(yTempRkck, ak2); // t + a2h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yTempRkck[i][j]=y_[i][j]+h*(b31*dydt_[i][j]+b32*ak2[i][j]);
  
  O_->derivs(yTempRkck, ak3); // t + a3h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yTempRkck[i][j]=y_[i][j]+h*(b41*dydt_[i][j]+b42*ak2[i][j]+b43*ak3[i][j]);
  
  O_->derivs(yTempRkck, ak4); // t + a4 * h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yTempRkck[i][j]=y_[i][j]+h*
				(b51*dydt_[i][j]+b52*ak2[i][j]+b53*ak3[i][j]+b54*ak4[i][j]);
  
  O_->derivs(yTempRkck, ak5); // t + a5 * h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yTempRkck[i][j]=y_[i][j]+h*
				(b61*dydt_[i][j]+b62*ak2[i][j]+b63*ak3[i][j]+b64*ak4[i][j]+
				 b65*ak5[i][j]);
  
  O_->derivs(yTempRkck, ak6); // t + a6 * h
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			yOut[i][j]=y_[i][j]+h*
				(c1*dydt_[i][j]+c3*ak3[i][j]+c4*ak4[i][j]+c6*ak6[i][j]);
  
	for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      yErr[i][j]=h*
        (dc1*dydt_[i][j]+dc3*ak3[i][j]+dc4*ak4[i][j]+
         dc5*ak5[i][j]+dc6*ak6[i][j]);
}

void RK5AdaptiveEquilibrium::
rkckMechanical(double h,
	       std::vector< std::vector<double> > &yOut,
	       std::vector< std::vector<double> > &yErr,
	       std::vector< std::vector<double> > &ak2,
	       std::vector< std::vector<double> > &ak3,
	       std::vector< std::vector<double> > &ak4,
	       std::vector< std::vector<double> > &ak5,
	       std::vector< std::vector<double> > &ak6,
	       std::vector< std::vector<double> > &yTempRkck,
	       std::vector<size_t> &mechEq,std::vector<size_t> &posVar ) 
{  
  //static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875;
  static double b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+b21*h*dydt_[i][j];
    }
  O_->derivsMechanical(yTempRkck,ak2,mechEq,posVar);//t+a2h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+h*(b31*dydt_[i][j]+b32*ak2[i][j]);
    }
  O_->derivsMechanical(yTempRkck,ak3,mechEq,posVar);//t+a3h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+h*(b41*dydt_[i][j]+b42*ak2[i][j]+b43*ak3[i][j]);
    }
  O_->derivsMechanical(yTempRkck,ak4,mechEq,posVar);//t+a4*h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+h*
        (b51*dydt_[i][j]+b52*ak2[i][j]+b53*ak3[i][j]+b54*ak4[i][j]);
    }
  O_->derivsMechanical(yTempRkck,ak5,mechEq,posVar);//t+a5*h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yTempRkck[i][j]=y_[i][j]+h*
        (b61*dydt_[i][j]+b62*ak2[i][j]+b63*ak3[i][j]+b64*ak4[i][j]+
         b65*ak5[i][j]);
    }  
  O_->derivsMechanical(yTempRkck,ak6,mechEq,posVar);//t+a6*h
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yOut[i][j]=y_[i][j]+h*
        (c1*dydt_[i][j]+c3*ak3[i][j]+c4*ak4[i][j]+c6*ak6[i][j]);
    }
  for(size_t i=0 ; i<N() ; i++ )
    for(size_t k=0 ; k<posVar.size() ; k++ ) {
      size_t j=posVar[k];
      yErr[i][j]=h*
        (dc1*dydt_[i][j]+dc3*ak3[i][j]+dc4*ak4[i][j]+
         dc5*ak5[i][j]+dc6*ak6[i][j]);
    }
}




//////////////////////////////////////////////////////////////////////
//                                                                  //
//               CLASS SOLVERRK5ADAPTIVETEMPLATE                    //
//                                                                  //
//////////////////////////////////////////////////////////////////////


RK5AdaptiveTemplate::RK5AdaptiveTemplate(Organism *O,std::ifstream &IN)
	:BaseSolver(O,IN)
{
	readParameterFile(IN);
}

RK5AdaptiveTemplate::~RK5AdaptiveTemplate()
{
	//std::cerr << "delete costFunction...";
	delete C_;
	//std::cerr << "Done!\n";
}

void RK5AdaptiveTemplate::readParameterFile(std::ifstream &IN, int verbose)
{
  verbose=1;
  if (verbose) {
    std::cerr << "RK5AdaptiveTemplate: Reading parameter file:" << std::endl;
  }
  std::string costType;
  IN >> costType;
  C_ = BaseCostFunction::createCostFunction(costType);
  
  IN >> startTime_;
  t_= startTime_;
  IN >> endTime_;
  
  IN >> printFlag_;
  IN >> numPrint_;
  
  IN >> h1_;
  IN >> eps_;
  
  if (verbose) {
    std::cerr << "Cost function: " << costType << std::endl
	      << "T_start/end: " << startTime_ << " " << endTime_ << std::endl
	      << "Print_flag: " << printFlag_ << " PrintNum: " << numPrint_ << std::endl
	      << "Initial/max step: " << h1_ << " eps: " << eps_ << std::endl;
  }
  
  //simulation template parameters
  size_t Mtmp;
  IN >> Mtmp;
  
  simulationIndex_.resize( Mtmp );
  
  if (Mtmp) {
    if (verbose)
      std::cerr << Mtmp << " variables read from template." << std::endl;
    for (size_t i=0 ; i<Mtmp ; ++i) {
      IN >> simulationIndex_[i];
      if (verbose)
	std::cerr << simulationIndex_[i] << " ";
    }
    IN >> inputTemplateFile_;
    if (verbose)
      std::cerr << std::endl << "Template file " << inputTemplateFile_ << std::endl; 
  }
  else
    if (verbose)
      std::cerr << "No simulation template file used." << std::endl; 
  
  //cost template parameters
  IN >> Mtmp;
  std::vector<size_t> costList;
  
  if( Mtmp ) {
    if (verbose)
      std::cerr << Mtmp << " variables used in cost calculation\n";
    costList.resize(Mtmp);
    for (size_t i=0 ; i<Mtmp ; i++ ) {
      IN >> costList[i];
      if (verbose)
	std::cerr << costList[i] << " ";
    }
    if (verbose)
      std::cerr << "\n";
    setCostList(costList);
    std::string costTemplateFile;	
    IN >> costTemplateFile;
    setCostTemplateFile(costTemplateFile);
    if (verbose)
      std::cerr << "Cost template file " << costTemplateFile << "\n"; 
    C_->initiateCostCalculation(startTime_, endTime_);
  }
  else
    if (verbose)
      std::cerr << "No cost template file used\n";
  //IN.close();
  if (verbose) {
    std::cerr << "RK5AdaptiveTemplate: Done reading parameter file." << std::endl << std::endl;
  }
}

void RK5AdaptiveTemplate::simulate(void)
{ 
	static bool flag = true;
	if(flag)
		setSimulationFlag();
	flag = false;
	
  // double tiny = 1e-30; // Using NR definition
  double tiny = 1e-9*eps_; // Caveat! Using new definition
	double  h, hNext, hDid;
  //Vectors to store division data from template file
  std::vector<int> divIndex;
  std::vector<double> divTime;
  size_t division=0;

  //Check if template reading is needed...
  size_t simCount=0;
  size_t templateFlag=1;
  for (size_t i=0 ; i<M() ; i++ ) {
    simCount += simulationFlag_[i];
  }
  if( simCount==M() ) 
    templateFlag = 0;
  //Check if the costList is well defined
  for (size_t i=0; i<C_->sizeOfCostList(); ++i) {
    if (C_->costList(i)>=M()) {
      std::cerr << "RK5AdaptiveTemplate::simulate(): "
		<< C_->costList(i) << " is out of scope\n";
      exit(-1);
    }
  }
  // Check that h1 and endTime - startTime are > 0
  // 
  if (h1_ > 0.0 && (endTime_ - startTime_) > 0.0)
    h = h1_;
  else {//either h or (endTime-startTime) <=0
    std::cerr << "solverRk5Adaptive::simulate() - "
	      << "Wrong time borders or time step for simulation. "
	      << "No simulation performed.\n";
    exit(-1);
  }

  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  // Introduce the neighborhood
  // 
  if (O_->numNeighborhood())
    O_->neighborhoodCreate(y_, t_);
  if (y_.size() && y_.size() != dydt_.size()) {
    dydt_.resize(y_.size(), y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
  //Create all vectors that will be needed here and by rkqs and rkck!
  //////////////////////////////////////////////////////////////////////
  //Used here
  std::vector< std::vector<double> > yScal( N() );
  //Used by rkqs
  std::vector< std::vector<double> > yTemp( N() ),yErr( N() );
  //used by rkck
  std::vector< std::vector<double> > ak2( N() ),ak3( N() ),ak4( N() ),
    ak5( N() ),ak6( N() ),yTempRkck( N() );
  //Resize each vector
  for (size_t i = 0; i < N(); i++) {
    yScal[i].resize(M());
    yErr[i].resize(M());
    yTemp[i].resize(M());
    ak2[i].resize(M());
    ak3[i].resize(M());
    ak4[i].resize(M());
    ak5[i].resize(M());
    ak6[i].resize(M());
    yTempRkck[i].resize(M());
  }
  assert( y_.size()==yScal.size() && y_.size()==yErr.size() );
  assert( y_.size()==yTemp.size() && y_.size()==ak2.size() );
  assert( y_.size()==ak3.size() && y_.size()==ak4.size() );
  assert( y_.size()==ak5.size() && y_.size()==ak6.size() );
  assert( y_.size()==yTempRkck.size() );
  
  
  // Initiate cost Calculation if applicable
  //////////////////////////////////////////////////////////////////////
  double endTimeLocal = endTime_;
  if( C_->sizeOfCostList() ) {
    //C_->initiateCostCalculation(startTime_, endTime_);
    C_->resetCostCalculation();
    if (C_->sizeOfCostTemplateTime() ) {
      endTimeLocal = C_->endTime() +tiny;  
    }
  }
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
  double printTime = endTimeLocal + tiny;
  double printDeltaTime = endTimeLocal + 2.0 * tiny;
  if (numPrint_ <= 0) //No printing
    printFlag_ = 0;
  else if (numPrint_ == 1) { // Print last point (default)
  }
  else if (numPrint_ == 2) { //Print first/last point
    printTime = startTime_ - tiny;
  } 
  else { //Print first/last points and spread the rest uniformly
    printTime = startTime_ - tiny;
    printDeltaTime = (endTimeLocal - startTime_) / ((double) (numPrint_ - 1));
  }
  
  // Go
  //////////////////////////////////////////////////////////////////////
  t_ = startTime_;
  // Initiate template from file if applicable
  //////////////////////////////////////////////////////////////////////
  std::ifstream IN;
  if( templateFlag ) {
    const char *tmp = inputTemplateFile_.c_str();
    IN.open( tmp );
    if( !IN ) {
      std::cerr << "Rk5::simulateTemplate() "
		<< "Cannot open templatefile " << inputTemplateFile_
		<< "\n\n\7";exit(-1);}
    //Read previous and future time points
    initiateValuesTemplate(IN,divIndex,divTime);
  }
  
  numOk_ = numBad_ = 0;
	for (unsigned int nstp = 0;; nstp++) {
		if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
		// Update the derivatives
		O_->derivsTemplate(y_,dydt_,simulationFlag_);

		// Calculate 'scaling' for error measure
		for (size_t i = 0; i < N(); i++)
			for (size_t j = 0; j < M(); j++) {
				yScal[i][j] = fabs(y_[i][j]) + fabs(dydt_[i][j] * h) + tiny;
				// Try this modification to see whether it fasten up things
				// if( yScal[i][j]<1.0 )
				// yScal[i][j] = 1.0;
			}

		// Print if applicable 
		if (printFlag_ && t_ >= printTime) {
			printTime += printDeltaTime;
			print();
		}

		// Check if step is larger than max allowed
		// max step end is min of endTime_ and printTime
		double tMin = endTimeLocal < printTime ? endTimeLocal : printTime;

		int maxStepFlag = 0;//1 if stepsize shortened
		double hTemp=h; //stepsize sugested before shortened

		if( templateFlag && futureTime()<tMin ) tMin = futureTime(); 
		if( templateFlag && division<divIndex.size() && divTime[division]<tMin ) 
			tMin=divTime[division];
		if( C_->sizeOfCostList() && C_->costTemplateIndex()<C_->sizeOfCostTemplateTime() 
				&& C_->costTemplateTime( C_->costTemplateIndex() )<tMin )
			tMin=C_->costTemplateTime( C_->costTemplateIndex() );
		if( t_+h>tMin ){
			h=tMin-t_;
			maxStepFlag = 1;
		}

		// Update
		rkqsTemplate(h, hDid, hNext, yScal, yTemp, yErr, ak2, ak3, ak4, ak5, ak6, 
				yTempRkck);
		if (hDid == h) ++numOk_; else ++numBad_;
		//Update other cols of y from template
		if( templateFlag )
			allValuesTemplate();//time is already updated

		// Make updates of rections that requires it
		O_->reactionUpdate(hDid, t_, y_);

		// Check for cell divisions, compartmental removal etc
		//O_->compartmentChange(y_, dydt_, t_);
		if( templateFlag ) {
			int tmpDivCount=0;
			while( division<divIndex.size() && t_>=divTime[division] ) {
				divideTemplate(divIndex[division]);
				division++;
				tmpDivCount++;
			}
			if( tmpDivCount ) {
				// Rescale all temporary vectors as well
				if (y_.size() != yScal.size()) {
					yScal.resize(N(), yScal[0]);
					yTemp.resize(N(), yTemp[0]);
					yErr.resize(N(), yErr[0]);
					ak2.resize(N(), ak2[0]);
					ak3.resize(N(), ak3[0]);
					ak4.resize(N(), ak4[0]);
					ak5.resize(N(), ak5[0]);
					ak6.resize(N(), ak6[0]);
					yTempRkck.resize(N(), yTempRkck[0]);
				}
			}
		}
		//Add cost from template comparison if applicable
		//////////////////////////////////////////////////////////////////////
		//size_t tmpCount=0;
		//if (!tmpCount) {
		//std::cerr << "Cost calculation variables:" << std::endl;
		//std::cerr << C_->sizeOfCostList() << " "
		//		<< C_->costTemplateIndex() << " "
		//		<< C_->sizeOfCostTemplateTime() << " "
		//		<< C_->costTemplateTime(C_->costTemplateIndex()) << std::endl;
		//++tmpCount;
		//}
		if( C_->sizeOfCostList() && C_->costTemplateIndex()< C_->sizeOfCostTemplateTime()
				&& t_>= C_->costTemplateTime(C_->costTemplateIndex() ) ) {
			//std::cerr << "Cost added at time " << t_ << std::endl;
			C_->addCost(y_,t_); 
		}

		// Check if new time point data needs to be read from file
		//////////////////////////////////////////////////////////////////////
		if( templateFlag && t_>=futureTime_ && t_<endTimeLocal ) {
			updateValuesTemplate(IN,endTimeLocal,divIndex,divTime);
			division=0;    
		}
		// Update neighborhood
		if (O_->numNeighborhood())
			O_->neighborhoodUpdate(y_, t_);

		// If the end t is passed return (print if applicable)
		if (t_ >= endTimeLocal) {
			if (printFlag_) {
				// Update the derivatives
				O_->derivs(y_, dydt_);
				print();
			}
			//std::cerr << "Simulation done.\n"; 
			return;
		}
		//Warn for small step sizes...
		//  if (fabs(hNext) <= hMin) {
		//        std::cerr << "Warning: Step size small (" << hNext
		//  		<< ") in rk5Adaptive::simulate at time " 
		//  		<< t << "\n"; /*exit(-1);*/ }
		if(!maxStepFlag)
			h=hNext;
		else
			h=hTemp;
		//Do not take larger steps than h1
		if (h > h1_)
			h = h1_;
	}
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
void RK5AdaptiveTemplate::rkqsTemplate(double hTry,double &hDid,double &hNext,
		std::vector< std::vector<double> > &yScal,
		std::vector< std::vector<double> > &yTemp,
		std::vector< std::vector<double> > &yErr,
		std::vector< std::vector<double> > &ak2,
		std::vector< std::vector<double> > &ak3,
		std::vector< std::vector<double> > &ak4,
		std::vector< std::vector<double> > &ak5,
		std::vector< std::vector<double> > &ak6,
		std::vector< std::vector<double> > &yTempRkck ) 
{  
	double errMax,h,hTemp,tNew,aux;
	h=hTry;
	double nanProtection = 10000.*eps_; // 2.*eps_; // must be larger than 1*eps_
	size_t numTopologyVariable = O_->numTopologyVariable();
	//rkckTemplate(h,yTemp,yErr,
	// ak2,ak3,ak4,ak5,ak6,yTempRkck);
	for (;;) {
		rkckTemplate(h,yTemp,yErr,
				ak2,ak3,ak4,ak5,ak6,yTempRkck);
		errMax=0.;
		for (size_t j=0 ; j<M() ; j++ )
			if( simulationFlag_[j] )
				for (size_t i=0 ; i<N() ; i++ ) {
					aux = fabs(yErr[i][j]/yScal[i][j]);
					if( aux>errMax )
						errMax = aux;
					else if (std::isnan(aux)){ // if aux==nan, ensure that we reduce step size.
						//std::cerr << "NAN forced stepsize decrease in rungeKutta.cc RK5AdaptiveTemplate::rkqs\n";
						if (errMax <  nanProtection)
							errMax = nanProtection;
					}
					if( yTemp[i][j] < 0. && j >= numTopologyVariable){ // Reduce step size if the proposed step resulted in a negative concentration.
						if (errMax < nanProtection)
							errMax = nanProtection;
					}
				}
		errMax /= eps_;
		if (errMax <= 1.0) break;
		//if (errMax <= 1.0 || h <= 0.01) break;
		hTemp=SAFETY*h*pow(errMax,PSHRNK);
		if( h >= 0. )
			h = hTemp > 0.1*h ? hTemp : 0.1*h;
		else
			h = hTemp > 0.1*h ? 0.1*h : hTemp;
		tNew=t_+h;
		if( tNew==t_ ) { std::cerr << "stepsize underflow in RK5AdaptiveTemplate::rkqsTemplate\n"
			<< "This may for example happen if you have a negative concentration.\n" ; exit(-1); }
	}
	if (errMax  <= 1.0) {  // why is this here? the only way to get here is by beaking the above loop. This cannot evaluate as false.
		if (errMax > ERRCON) hNext=SAFETY*h*pow(errMax,PGROW);
		else hNext=5.0*h;
	}
	else
		hNext = h;
	//std::cerr << h << " " << hTry << " " << hNext << " " << hDid << " " << "\n";
	t_ += (hDid=h);
	for (size_t j=0 ; j<M() ; j++ )
		if( simulationFlag_[j] )
			for (size_t i=0 ; i<N() ; i++ )
				y_[i][j]=yTemp[i][j];

}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

void RK5AdaptiveTemplate::rkckTemplate(double h,
		       std::vector< std::vector<double> > &yOut,
		       std::vector< std::vector<double> > &yErr,
		       std::vector< std::vector<double> > &ak2,
		       std::vector< std::vector<double> > &ak3,
		       std::vector< std::vector<double> > &ak4,
		       std::vector< std::vector<double> > &ak5,
		       std::vector< std::vector<double> > &ak6,
		       std::vector< std::vector<double> > &yTempRkck ) 
{  
    //static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875;
    static double b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  
  for (size_t i=0 ; i<N() ; i++ )
    for (size_t j=0 ; j<M() ; j++ )
      yTempRkck[i][j]=y_[i][j]+b21*h*dydt_[i][j];
  
  O_->derivsTemplate(yTempRkck,ak2,simulationFlag_);//t+a2h
  for (size_t i=0 ; i<N() ; i++ )
    for (size_t j=0 ; j<M() ; j++ )
      yTempRkck[i][j]=y_[i][j]+h*(b31*dydt_[i][j]+b32*ak2[i][j]);
  
  O_->derivsTemplate(yTempRkck,ak3,simulationFlag_);//t+a3h
  for (size_t i=0 ; i<N() ; i++ )
    for (size_t j=0 ; j<M() ; j++ )
      yTempRkck[i][j]=y_[i][j]+h*(b41*dydt_[i][j]+b42*ak2[i][j]+b43*ak3[i][j]);
  
  O_->derivsTemplate(yTempRkck,ak4,simulationFlag_);//t+a4*h
  for (size_t i=0 ; i<N() ; i++ )
    for (size_t j=0 ; j<M() ; j++ )
      yTempRkck[i][j]=y_[i][j]+h*
        (b51*dydt_[i][j]+b52*ak2[i][j]+b53*ak3[i][j]+b54*ak4[i][j]);
  
  O_->derivsTemplate(yTempRkck,ak5,simulationFlag_);//t+a5*h
  for (size_t i=0 ; i<N() ; i++ )
    for (size_t j=0 ; j<M() ; j++ )
      yTempRkck[i][j]=y_[i][j]+h*
        (b61*dydt_[i][j]+b62*ak2[i][j]+b63*ak3[i][j]+b64*ak4[i][j]+
         b65*ak5[i][j]);
  
  O_->derivsTemplate(yTempRkck,ak6,simulationFlag_);//t+a6*h
  for (size_t j=0 ; j<M() ; j++ )
    if( simulationFlag_[j] )
      for (size_t i=0 ; i<N() ; i++ )
	yOut[i][j]=y_[i][j]+h*
	  (c1*dydt_[i][j]+c3*ak3[i][j]+c4*ak4[i][j]+c6*ak6[i][j]);
  
  for (size_t j=0 ; j<M() ; j++ )
    if( simulationFlag_[j] )
      for (size_t i=0 ; i<N() ; i++ )
	yErr[i][j]=h*
	  (dc1*dydt_[i][j]+dc3*ak3[i][j]+dc4*ak4[i][j]+
	   dc5*ak5[i][j]+dc6*ak6[i][j]);
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
//                          CLASS SOLVERRK4                          //
//                                                                   //
///////////////////////////////////////////////////////////////////////

RK4::RK4(Organism *O,std::ifstream &IN)
	:BaseSolver(O,IN)
{
	readParameterFile(IN);
}

void RK4::readParameterFile(std::ifstream &IN)
{
  //Read in the needed parameters
  IN >> startTime_;
  t_=startTime_;
  IN >> endTime_;
  
  IN >> printFlag_;
  IN >> numPrint_;
  
  IN >> h_;
}

void RK4::simulate() 
{
  // Check that h1 and endTime-startTime are > 0
  //
  if( !(h_>0. && (endTime_-startTime_)>0.) ) {
    std::cerr << "Rk4::simulate() Wrong time borders or time step for "
	      << "simulation. No simulation performed.\n";
    return;
  }
	std::cerr << "Simulating using fourth-order Runge-Kutta\n";

  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  // Introduce the neighborhood
  //
  if( O_->numNeighborhood() )
    O_->neighborhoodCreate(y_,t_);
  if( y_.size() && y_.size() != dydt_.size() ) {
    dydt_.resize( y_.size(),y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
	
  // Create all vectors that will be needed by rk4()!
  //
  std::vector< std::vector<double> > yt( N() ),dyt( N() ),dym( N() );
  //Resize each vector
  for( size_t i=0 ; i<N() ; i++ ) {
    yt[i].resize(M());
    dyt[i].resize( M() );
    dym[i].resize( M() );
  }
  assert( y_.size()==yt.size() && y_.size()==dyt.size() &&
	  y_.size()==dym.size() );
  
  // Initiate print times
  //
  double tiny = 1e-10;
  double printTime=endTime_+tiny;
  double printDeltaTime=endTime_+2.*tiny;
  if( numPrint_<=0 )//No printing
    printFlag_=0;
  else if( numPrint_==1 ) {//Print last point (default)
  }
  else if( numPrint_==2 ) {//Print first/last point
    printTime=startTime_-tiny;
  } 
  else {//Print first/last points and spread the rest uniformly
    printTime=startTime_-tiny;
    printDeltaTime=(endTime_-startTime_)/double(numPrint_-1);
  }
  
  // Go
  //
  t_=startTime_;
  numOk_ = numBad_ = 0;
  while( t_<endTime_ ) {
		if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
   //Update the derivatives
    O_->derivs(y_,dydt_);
    
    //Print if applicable 
    if( printFlag_ && t_ >= printTime ) {
      printTime += printDeltaTime;
      print();
    }

    //Check if step is larger than max allowed
    //max step end is min of endTime_ and printTime
    //double tMin= endTime_<printTime ? endTime_ : printTime;
    //if( t_+h>tMin ) h=tMin-t_;

    //Update
    rk4(yt,dyt,dym);
    numOk_++;
    
		// Make updates of rections that requires it
		O_->reactionUpdate(h_, t_, y_);

    // Check for cell divisions, compartmental removal etc
        //DISABLE 042521
    //O_->compartmentChange(y_,dydt_,t_);
    
    //Rescale all temporary vectors as well
    if(y_.size() != yt.size() ) {
      yt.resize( N(), yt[0] );
      dyt.resize( N(), dyt[0] );
      dym.resize( N(), dym[0] );
    }
    
    //Update neighborhood
    if( O_->numNeighborhood() )
      O_->neighborhoodUpdate(y_,t_);
    
    //update time variable
    if( (t_+h_)==t_ ) {
      std::cerr << "Rk4::simulate() Step size too small.";
      exit(-1);
    }
    t_ += h_;
  }
  if( printFlag_ ) {
    //Update the derivatives
    O_->derivs(y_,dydt_);
    print();
  }
  std::cerr << "Simulation done.\n"; 
  return;
}

void RK4::simulateWithConvergingMechanics(double cEps,double cTmax,
					  std::vector<size_t> &mechEq,
					  std::vector<size_t> &mechPos) 
{
  //Check that h1 and endTime-startTime are > 0
  //////////////////////////////////////////////////////////////////////
  if( !(h_>0. && (endTime_-startTime_)>0.) ) {
    std::cerr << "Rk4::simulate() Wrong time borders or time step for "
	      << "simulation. No simulation performed.\n";
    return;
  }
  
  // Initiate reactions for those where it is applicable
  //
  O_->reactionInitiate(startTime_,y_);

  //Introduce the neighborhood
  //////////////////////////////////////////////////////////////////////
  if( O_->numNeighborhood() )
    O_->neighborhoodCreate(y_,t_);
  if( y_.size() && y_.size() != dydt_.size() ) {
    dydt_.resize( y_.size(),y_[0]);
  }
  assert( y_.size() == O_->numCompartment() && y_.size()==dydt_.size());
  
  //Create all vectors that will be needed by rk4()!
  //////////////////////////////////////////////////////////////////////
  std::vector< std::vector<double> > yt( N() ),dyt( N() ),dym( N() );
  //Resize each vector
  for( size_t i=0 ; i<N() ; i++ ) {
    yt[i].resize(M());
    dyt[i].resize( M() );
    dym[i].resize( M() );
  }
  assert( y_.size()==yt.size() && y_.size()==dyt.size() &&
	  y_.size()==dym.size() );
  
  // Initiate print times
  //////////////////////////////////////////////////////////////////////
  double tiny = 1e-10;
  double printTime=endTime_+tiny;
  double printDeltaTime=endTime_+2.*tiny;
  if( numPrint_<=0 )//No printing
    printFlag_=0;
  else if( numPrint_==1 ) {//Print last point (default)
  }
  else if( numPrint_==2 ) {//Print first/last point
    printTime=startTime_-tiny;
  } 
  else {//Print first/last points and spread the rest uniformly
    printTime=startTime_-tiny;
    printDeltaTime=(endTime_-startTime_)/double(numPrint_-1);
  }
  
  // Go
  //////////////////////////////////////////////////////////////////////
  t_=startTime_;
  numOk_ = numBad_= 0;
  while( t_<endTime_ ) {
    if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
    //Update the derivatives
    O_->derivs(y_,dydt_);
    
    //Print if applicable 
    if( printFlag_ && t_ >= printTime ) {
      printTime += printDeltaTime;
      print();
    }

    //Check if step is larger than max allowed
    //max step end is min of endTime_ and printTime
    //double tMin= endTime_<printTime ? endTime_ : printTime;
    //if( t_+h>tMin ) h=tMin-t_;

    //Update
    rk4(yt,dyt,dym);
    numOk_++;

    //Converge the mechanics
    updateMechanicsUntilConvergence(cEps,cTmax,yt,dyt,dym,mechEq,mechPos);

		// Make updates of rections that requires it
		O_->reactionUpdate(h_, t_, y_);

    // Check for cell divisions, compartmental removal etc
        //DISABLED 042521
    //O_->compartmentChange(y_,dydt_,t_);
    
    //Rescale all temporary vectors as well
    if(y_.size() != yt.size() ) {
      yt.resize( N(), yt[0] );
      dyt.resize( N(), dyt[0] );
      dym.resize( N(), dym[0] );
    }
    
    //Update neighborhood
    if( O_->numNeighborhood() )
      O_->neighborhoodUpdate(y_,t_);
    
    //update time variable
    if( (t_+h_)==t_ ) {
      std::cerr << "Rk4::simulate() Step size too small.";
      exit(-1);
    }
    t_ += h_;
  }
  if( printFlag_ ) {
    //Update the derivatives
    O_->derivs(y_,dydt_);
    print();
  }
  std::cerr << "Simulation done.\n"; 
  return;
}

void RK4::updateMechanicsUntilConvergence(double cEps,double cTmax,
					  std::vector< std::vector<double> > 
					  &yt,
					  std::vector< std::vector<double> > 
					  &dyt,
					  std::vector< std::vector<double> > 
					  &dym,
					  std::vector<size_t> &mechEq,
					  std::vector<size_t> &mechPos ) 
{ 
  // Go
  //////////////////////////////////////////////////////////////////////
  double t = 0.0, maxDydt=0.0;
  while( t<cTmax ) {
    if (debugFlag()) {
			yCopy_[debugCount()] = y_;
		} 
    //Update the derivatives
    O_->derivsMechanical(y_,dydt_,mechEq,mechPos);
    
    //Update
    rk4Mechanical(yt,dyt,dym,mechEq,mechPos);
    
		// Make updates of rections that requires it
		O_->reactionUpdate(h_, t_, y_);

    // Check for cell divisions, compartmental removal etc
        //DISABLED 042521
    //O_->compartmentChange(y_,dydt_,t_);
    
    //Rescale all temporary vectors as well
    if(y_.size() != yt.size() ) {
      yt.resize( N(), yt[0] );
      dyt.resize( N(), dyt[0] );
      dym.resize( N(), dym[0] );
    }
    
    //Update neighborhood
    if( O_->numNeighborhood() )
      O_->neighborhoodUpdate(y_,t_);
    
    //Check if derivatives have converged
    maxDydt=0.0;
    for( size_t i=0 ; i<N() ; i++ )
      for( size_t k=0 ; k<mechPos.size() ; k++ ) {
	size_t j=mechPos[k];
	if( std::fabs(dydt_[i][j])>maxDydt )
	  maxDydt = std::fabs(dydt_[i][j]); 
      }
    if( maxDydt<cEps ) {
      //std::cerr << "Rk4::updateMechanicsUntilConvergance "
      //	<< "t=" << t << " maxDydt=" << maxDydt << "\n";
      return;
    }
    
    
    t += h_;
  }
  std::cerr << "Rk4::updateMechanicsUntilConvergance "
	    << "Warning max time " << cTmax << " reached with"
	    << " Max dy/dt " << maxDydt << "." << std::endl; 
  return;
}

void RK4::rk4(std::vector< std::vector<double> > &yt,
	      std::vector< std::vector<double> > &dyt,
	      std::vector< std::vector<double> > &dym ) 
{  
  double hh=0.5*h_;
  double h6=h_/6.0;
	size_t n = N();
	size_t m = M();
  for(size_t i=0 ; i<n ; ++i )
    for( size_t j=0 ; j<m ; ++j )
      yt[i][j]=y_[i][j]+hh*dydt_[i][j];
  
  O_->derivs(yt,dyt);//second step
  for(size_t i=0 ; i<n ; ++i )
    for( size_t j=0 ; j<m ; ++j )
      yt[i][j]=y_[i][j]+hh*dyt[i][j];
  
  O_->derivs(yt,dym);//third step
  for(size_t i=0 ; i<n ; ++i )
    for( size_t j=0 ; j<m ; ++j ) {
      yt[i][j]=y_[i][j]+h_*dym[i][j];
      dym[i][j] += dyt[i][j];
    }
  
  O_->derivs(yt,dyt);//fourth step
  for(size_t i=0 ; i<n ; ++i )
    for( size_t j=0 ; j<m ; ++j )
      y_[i][j]=y_[i][j]+h6*(dydt_[i][j]+dyt[i][j]+2.0*dym[i][j]);
}

void RK4::rk4Mechanical(std::vector< std::vector<double> > &yt,
			std::vector< std::vector<double> > &dyt,
			std::vector< std::vector<double> > &dym,
			std::vector<size_t> &mechEq,
			std::vector<size_t> &mechPos ) 
{  
  double hh=0.5*h_;
  double h6=h_/6.0;

  for(size_t i=0 ; i<N() ; i++ )
    for( size_t k=0 ; k<mechPos.size() ; k++ ) {
      size_t j=mechPos[k];
      yt[i][j]=y_[i][j]+hh*dydt_[i][j];
    }
  
  O_->derivsMechanical(yt,dyt,mechEq,mechPos);//second step
  for(size_t i=0 ; i<N() ; i++ )
    for( size_t k=0 ; k<mechPos.size() ; k++ ) {
      size_t j=mechPos[k];
      yt[i][j]=y_[i][j]+hh*dyt[i][j];
    }
  
  O_->derivsMechanical(yt,dym,mechEq,mechPos);//third step
  for(size_t i=0 ; i<N() ; i++ )
    for( size_t k=0 ; k<mechPos.size() ; k++ ) {
      size_t j=mechPos[k];
      yt[i][j]=y_[i][j]+h_*dym[i][j];
      dym[i][j] += dyt[i][j];
    }
  
  O_->derivsMechanical(yt,dyt,mechEq,mechPos);//fourth step
  for(size_t i=0 ; i<N() ; i++ )
    for( size_t k=0 ; k<mechPos.size() ; k++ ) {
      size_t j=mechPos[k];
      y_[i][j]=y_[i][j]+h6*(dydt_[i][j]+dyt[i][j]+2.0*dym[i][j]);
    }
}

//!Finds the maximal |dydt|/|y| for the system
double RK4::maxDerivative() {  
  double max=0.0,val;
  for( size_t i=0 ; i<N() ; i++ ) {
    for( size_t j=0 ; j<M() ; j++ ) {
      if( y_[i][j]!=0 )
	val = fabs( dydt_[i][j]/y_[i][j] );
      else
	val = fabs( dydt_[i][j] );
      
      if( val>max ) 
	max = val;
    }
  }
  return max;
}

