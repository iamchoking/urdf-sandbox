# URDF Sandbox

Collection of small RaiSim demos for experimenting with RAIPAL and other URDF configurations.

## Dependencies
- [RaiSim](https://raisim.com/) (match the version used inside this workspace and ensure `raisim::raisim` is discoverable through `CMAKE_PREFIX_PATH`).
- [`raipal_kinematics`](https://github.com/railabatkaist/raipal_kinematics) – provides the crossed four-bar helpers shared with the production RAIPAL stack. Install it system-wide (e.g. `git clone ... && ./install.sh`) so `find_package(raipal_kinematics CONFIG REQUIRED)` succeeds.

Once both dependencies are installed, configure and build as usual:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/raisim
cmake --build build
```

`raipal_4-bar.cpp` will automatically link against `raipal_kinematics` and reuse its FK utilities instead of keeping duplicate polynomial code in this repository.
