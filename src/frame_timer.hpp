#ifndef FRAME_TIMER_HPP
#define FRAME_TIMER_HPP

#include <chrono>
#include <thread>
#include <iostream>

class FrameTimer {
public:
  FrameTimer(double timestep): timestep_(std::chrono::microseconds(static_cast<long>(timestep * 1e6))), frame_count_(0), started_(false), enabled_(true) {}

  void enable(double timestep)  {setTimestep(timestep); enable();}
  void enable()  {enabled_ = true ;started_ = false;}
  void disable() {enabled_ = false;started_ = false;}

  /// @brief starts the frame timer (DO NOT call this manually if using tick())
  void start() {
    if(!enabled_ || started_) return;
    loop_start_ = std::chrono::high_resolution_clock::now();
    elapsed_ = std::chrono::microseconds(0);
    nominal_ = std::chrono::microseconds(0);
    frame_count_ = 0;
    started_ = true;
  }

  /// @brief putting this method in the beginning of a for loop ensures a consistent framerate.
  void tick(){
    if(!enabled_) return;
    if(!started_) {start(); return;}
    frame_count_++;
    nominal_ += timestep_;

    elapsed_ = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - loop_start_);
    delta_ = nominal_ - elapsed_;

    if(delta_.count() > 0){
      std::this_thread::sleep_for(delta_);
    }
    else if (delta_.count() < -1e5){
      if(frame_count_ % 1000 == 0){
        std::cout << "[FrameTimer] Frame #" << frame_count_ << " is " << -delta_.count()/1000.0 << " ms slower than nominal." <<std::endl;
      }
    }
  }

  /// @brief manually enforce the next frame (for non-constant timestep scenarios)
  inline void tick(double timestep){setTimestep(timestep);tick();}

  /// @brief Ends frame timer and returns the total elapsed time in seconds.
  /// @return Delta time: (Positive: frames were faster than nominal)
  double end(){
    if(!enabled_) return 0.0;
    tick();
    started_ = false;
    elapsed_ = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - loop_start_);
    // std::cout << "[FrameTimer] Finished.\n\tTimestep: " << timestep_ << " s.\n\tTotal elapsed: " << elapsed_.count()/1e6 << " s." <<std::endl;
    return elapsed_.count() / 1.0e6;
  }

  /// @brief (re)sets nominal timestep
  /// @param timestep timestep in seconds
  inline void setTimestep(double timestep){
    timestep_ = std::chrono::microseconds(static_cast<long>(timestep * 1e6));
  }

private:

  inline void rest(){

  }

  bool enabled_,started_;
  std::chrono::microseconds timestep_;
  size_t frame_count_;
  std::chrono::microseconds elapsed_;
  std::chrono::microseconds nominal_;
  std::chrono::microseconds delta_;
  std::chrono::high_resolution_clock::time_point loop_start_;
};

#endif // FRAME_TIMER_HPP