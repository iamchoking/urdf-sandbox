#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"
#include <chrono>

#include "frame_timer.hpp"
#include "table_printer.hpp"

// double DT_MIN = 0.001;
double DT_MIN  = 0.0010;
double DT_MAX  = 0.0025;
double DT_STEP = 0.0001;

double CUTOFF_MIN  =  5.0;
double CUTOFF_MAX  = 50.0;
double CUTOFF_STEP =  1.0;

double SAMPLE_TIME = 2.5;
size_t NUM_SAMPLES = 20;

std::string URDF_TYPE = "_stub-0";

bool RANDOM_DT_STEP = true; // randomize dt every step

bool VISUALIZE = true;

// randomizers
std::mt19937 gen_;
std::normal_distribution<double> normDist_{0., 1.};
std::uniform_real_distribution<double> uniDist_{0., 1.};
std::exponential_distribution<double> simDtExpDist_{400.};
std::exponential_distribution<double> conDtExpDist_{100.};

/// @brief Crossed 4-bar linkage forward kinematics calculation. Calulates index 4/5 of gc based on index 3.
/// @param gc 9D generalized coordinate vector (with index 4/5 as 0)
void cfbFK(Eigen::VectorXd &gc){
  static double coeffGC5[17] = {7.487348861572753, -88.17268027058002, 470.7676589474978, -1507.949571146181, 3232.924175366478, -4903.587201465159, 5428.248858581333, -4464.681332539945, 2751.745273804874, -1265.589676288956, 418.4687974259759, -86.36058832217348, 6.540591054554072, -1.820824957906633, 2.514881046035373, 1.022698355939166, 0.000002629580961635315};
  static double coeffGC4[17] = {-4.599352352699746, 54.74231689702444, -295.2311950384566, 953.4731854791298, -2052.475077842654, 3099.83081925634, -3363.608943070999, 2636.542385991094, -1477.728449862623, 579.0562233124353, -157.9024440428828, 35.05094543077983, -7.773109367008759, -0.656408504334005, 0.798469056970581, 2.098468226080274, -0.000001063969759856979};
  static double l0 = 112.95, l1 = 60.0, l2 = 45.0, l3 = 85.47477864;
  static double a0 = -0.02249827895, b0 = 0.99988266, c0 = 2.722713633;

  gc[5] = 0;
  gc[4] = 0;
  double gc3_pow = 1.0;
  for(int i=16; i>=0; i--){
    gc[5] += coeffGC5[i] * gc3_pow;
    gc[4] += coeffGC4[i] * gc3_pow;
    gc3_pow *= gc[3];
  }

}

void sampleTarget(Eigen::VectorXd &pTarget, const std::vector<raisim::Vec<2UL>> jointLimits){
    for (int i = 0; i < jointLimits.size(); i++) {
      pTarget[i] = jointLimits[i][0] + uniDist_(gen_) * (jointLimits[i][1] - jointLimits[i][0]);
    }
    cfbFK(pTarget);
}

int main(int argc, char* argv[]) {
  gen_.seed(std::chrono::system_clock::now().time_since_epoch().count());
  // create raisim world
  raisim::World world; // physics world
  std::unique_ptr<raisim::RaisimServer> server; // visualization server
  
  std::string urdf_prefix = std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal" + URDF_TYPE;
  
  auto raipal_FT = world.addArticulatedSystem(urdf_prefix + "_R.urdf");
  raipal_FT->setName("raipal_FT");
  raipal_FT->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
  auto raipal_PD = world.addArticulatedSystem(urdf_prefix + "_L.urdf");
  raipal_PD->setName("raipal_PD");
  raipal_PD->setControlMode(raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE);

  raisim::ArticulatedSystemVisual *raipalTarget_FT, *raipalTarget_PD;
  raisim::Visuals *indicator;
  
  if (VISUALIZE){
    server = std::make_unique<raisim::RaisimServer>(&world);
    // auto indicator = server.addVisualSphere("indicator", 0.05, 1.0,1.0,0.0,0.8,"",true);
    indicator = server -> addVisualSphere("indicator", 0.05, 0.0,0.0,0.0,0.8,"",true);
    indicator->setPosition(Eigen::Vector3d(0.0, 0.0, 0.8));
    raipalTarget_FT = server->addVisualArticulatedSystem("Target_FT", urdf_prefix + "_R.urdf", 1.0,0.0,0.0,0.3);
    raipalTarget_PD = server->addVisualArticulatedSystem("Target_PD", urdf_prefix + "_L.urdf", 1.0,0.0,0.0,0.3);
  }

  // unpowered joint indices: 4/5

  // raipal_FT -> setComputeInverseDynamics(true);

  int gcDim, gvDim, nJoints = 7;   // unpowered joint indices: 4/5
  gcDim = raipal_FT->getGeneralizedCoordinateDim();
  gvDim = raipal_FT->getDOF();
  // std::cout << "gcDim: " << gcDim << ", gvDim: " << gvDim << std::endl;

  Eigen::VectorXd gc(gcDim), gv(gvDim);
  gc.setZero();gv.setZero();

  // low-pass filter implementation
  double lpf_alpha;
  Eigen::VectorXd gcRaw(gcDim), gvRaw(gvDim);

  /// initial (nominal) pose
  Eigen::VectorXd gc_init, gv_init, pTarget, dTarget;
  gc_init.setZero(gcDim) ; gv_init.setZero(gvDim);
  pTarget.setZero(gcDim); dTarget.setZero(gvDim);

  /// set pd gains
  Eigen::VectorXd jointPgain(gvDim), jointDgain(gvDim);
  jointPgain.setOnes();
  jointPgain *= 200;
  // jointPgain[1] = 1000.0; // adduction needs stronger gravity compensation
  jointDgain.setOnes();
  jointDgain *= 20.0;

  // unpowered joint indices: 4/5
  jointPgain[4] = 0.0;
  jointPgain[5] = 0.0;
  jointDgain[4] = 0.0;
  jointDgain[5] = 0.0;

  // jointDgain.tail(nJoints).setConstant(1);

  raipal_PD->setPdGains(jointPgain, jointDgain);

  raipal_FT->setState(gc_init, gv_init);
  raipal_PD->setState(gc_init, gv_init);
  raipal_FT->setGeneralizedForce(Eigen::VectorXd::Zero(gvDim));
  raipal_PD->setGeneralizedForce(Eigen::VectorXd::Zero(gvDim));

  FrameTimer ft(1.0);

  int num_dt     = int(((DT_MAX - DT_MIN) / DT_STEP) + 1);
  int num_cutoff = int(((CUTOFF_MAX - CUTOFF_MIN) / CUTOFF_STEP) + 1);
  std::cout << "Starting FT vs PD comparison test: ETA: ~ " << int((num_dt * num_cutoff * NUM_SAMPLES * SAMPLE_TIME) / 60.0) << " min " << std::endl;

  if(VISUALIZE){
    server->launchServer();
    server->focusOn(raipal_FT);
    for (int sec=3; sec>0; sec--){
      ft.tick();
      std::cout << "Starting in [" << sec << "]..." << std::endl;
    }
    ft.end();

    std::cout << "START!" << std::endl;
  }
  else{
    ft.disable();
  }

  bprinter::TablePrinter tpFailure(&std::cout);
  tpFailure.AddColumn("#",5);
  tpFailure.AddColumn("P Target",10);
  tpFailure.AddColumn("[FT] gc",10);
  tpFailure.AddColumn("[PD] gc",10);
  tpFailure.AddColumn("D Target",10);
  tpFailure.AddColumn("[FT] gv",10);
  tpFailure.AddColumn("[PD] gv",10);
  tpFailure.AddColumn("[FT] Force",10);
  tpFailure.AddColumn("[PD] Force",10);

  bprinter::TablePrinter tpResult(&std::cout);
  tpResult.AddColumn( "#(/" + std::to_string(num_dt * num_cutoff) + ")", 9);
  tpResult.AddColumn( "DT_Max (s)", 10);
  tpResult.AddColumn( "Cutoff(Hz)", 10);
  tpResult.AddColumn( "#Fail(/" + std::to_string(NUM_SAMPLES) + ")", 10);

  double dt = DT_MIN;
  double dtRange_min = DT_MIN;
  double dtRange_max = DT_MAX;

  double cutoff_freq = CUTOFF_MIN;

  double elapsed_time = 0.0, nominal_time = 0.0;
  Eigen::VectorXd gc_FT, gv_FT, gc_PD, gv_PD;

  size_t numFail = 0;

  size_t testCount = 0;
  tpResult.PrintHeader();
  // for(int dtIdx = 0; dtIdx < num_dt; dtIdx++){
    // dtRange_max = DT_MIN + (DT_MAX - DT_MIN) * (dtIdx) /(num_dt-1);
  dtRange_max = DT_MIN;
  while (dtRange_max <= DT_MAX + 1e-8){
    dt = dtRange_max;
    // for(int cutoffIdx = 0; cutoffIdx < num_cutoff; cutoffIdx++){
      //   cutoff_freq = CUTOFF_MIN + (CUTOFF_MAX - CUTOFF_MIN) * (cutoffIdx) / (num_cutoff - 1);
    cutoff_freq = CUTOFF_MIN;
    while(cutoff_freq <= CUTOFF_MAX + 1e-8){
      // std::cout << "[Test " << dtIdx * NUM_CUTOFF + cutoffIdx + 1 << "/" << NUM_DT * NUM_CUTOFF << "]: DT Range [" << dtRange_min << " ~ " << dtRange_max << "],  Cutoff: " << cutoff_freq << " Hz ";
      // std::flush(std::cout);
      
      for(int sampleIdx = 0; sampleIdx < NUM_SAMPLES; sampleIdx++){
        // sample a target
        sampleTarget(pTarget, raipal_FT->getJointLimits());

        if(VISUALIZE){
          raipalTarget_FT->setGeneralizedCoordinate(pTarget);
          raipalTarget_PD->setGeneralizedCoordinate(pTarget);
        } 

        world.setTimeStep(dt);

        // PD case: just set Pd target
        raipal_PD->setPdTarget(pTarget,dTarget);

        bool fail = false;
        bool declared = false;

        ft.reset();
        
        if(VISUALIZE){
          indicator->setColor(1.0,1.0,0.0,0.8); // make indicator yellow on start
        }
        world.setWorldTime(0.0);
        while (world.getWorldTime() < SAMPLE_TIME) {
          if(VISUALIZE){ft.tick(dt);}
          // std::cout << world.getWorldTime() << std::endl;
          // FT case
          // apply low-pass filter
          raipal_FT->getState(gcRaw, gvRaw);
          lpf_alpha = 2.0 * M_PI * cutoff_freq * dt / (2.0 * M_PI * cutoff_freq * dt + 1.0);
          gc = lpf_alpha * gcRaw + (1.0 - lpf_alpha) * gc;
          gv = lpf_alpha * gvRaw + (1.0 - lpf_alpha) * gv;

          // calculate feedforward torque
          Eigen::VectorXd tau(gvDim);
          tau.setZero();
          // simple pd control to compute feedforward torque
          for(int i=0; i<gcDim; i++){
            double p_err = pTarget[i] - gc[i];
            double d_err = dTarget[i] - gv[i];
            tau[i] = jointPgain[i] * p_err + jointDgain[i] * d_err;
          }
          raipal_FT->setGeneralizedForce(tau);

          if (server){server->lockVisualizationServerMutex();}
          world.integrate();
          if (server){server->unlockVisualizationServerMutex();}
      
          raipal_FT->getState(gc_FT, gv_FT);
          raipal_PD->getState(gc_PD, gv_PD);

          // check for anomalies
          (gc_FT - gc_PD).cwiseAbs().maxCoeff() > 1e-3 ? fail = true : fail = fail;
          (gv_FT - gv_PD).cwiseAbs().maxCoeff() > 1e-2 ? fail = true : fail = fail;
          if(ft.getElapsedTime() < 1.5){fail = false;} // ignore initial transient period
          else if(!fail && VISUALIZE){
            indicator->setColor(0.0,1.0,0.0,0.8); // make indicator green on success
          }

          if(fail && !declared){
            numFail++;
            declared = true;

            if(VISUALIZE){
              indicator->setColor(1.0,0.0,0.0,0.8); // make indicator red on fail
            }
            else{
              // no need to let it run further in non-visualization mode
              break;
            }

            // Eigen::VectorXd force_FT = raipal_FT->getGeneralizedForce().e();
            // Eigen::VectorXd force_PD = raipal_PD->getGeneralizedForce().e();
            // std::cout << "[FAIL] Time: " << ft.getElapsedTime() << " / dt: " << world.getTimeStep() << " / alpha: " << lpf_alpha << std::endl;

            // tpFailure.PrintHeader();
            // for (int i=0; i<gcDim; i++){
            //   tpFailure  << i
            //       << pTarget[i]
            //       << gcRaw[i]
            //       << gc_PD[i]
            //       << dTarget[i]
            //       << gvRaw[i]
            //       << gv_PD[i]
            //       << force_FT[i]
            //       << force_PD[i];
            // }
            // tpFailure.PrintFooter();

          }

          if(RANDOM_DT_STEP){
            dt = dtRange_min + (dtRange_max - dtRange_min) * uniDist_(gen_);
            world.setTimeStep(dt);
            // ft.step(world.getTimeStep());
          }
        }
        // fail ? std::cout << "x" : std::cout << "o";
        // std::flush(std::cout);
        // if(!fail){std::cout << "[PASS] Sample #" << sampleIdx << " passed." << std::endl;}
        nominal_time += world.getWorldTime();
        elapsed_time += ft.end();
      }
      testCount++;
      // std::cout <<  " -> " << numFail << "/" << NUM_SAMPLES << " failed." << std::endl;
      
      tpResult 
        << testCount
        << dtRange_max
        << cutoff_freq
        << numFail;

      // no longer need to try further failure cases (higher cutoff freq)
      if(numFail == NUM_SAMPLES){break;}        
      numFail = 0;
      cutoff_freq += CUTOFF_STEP;
    }
    tpResult.PrintFooter();
    numFail = 0;

    dtRange_max += DT_STEP;
  }
  tpResult.PrintFooter();

  std::cout << "Nominal Time: " << nominal_time << " seconds" << std::endl;
  std::cout << "Elapsed Time: " << elapsed_time << " seconds" << std::endl;

  if(VISUALIZE){server->killServer();}

  std::cout<<"SIMULATION COMPLETE"<<std::endl;
  
  return 0;
}
