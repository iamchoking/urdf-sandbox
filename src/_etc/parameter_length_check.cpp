#include "raisim/World.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

#define STRINGIFY_DETAIL(x) #x
#define STRINGIFY(x) STRINGIFY_DETAIL(x)

namespace {

struct ModuleCase {
  std::string target;
  fs::path seedPath;
  fs::path packageRoot;
  std::vector<std::string> modules;
};

struct GeneratedModuleConfig {
  fs::path configPath;
  std::string modulesForWorld;
};

std::string trim(const std::string &text) {
  const auto first = std::find_if_not(text.begin(), text.end(), [](unsigned char c) { return std::isspace(c); });
  const auto last = std::find_if_not(text.rbegin(), text.rend(), [](unsigned char c) { return std::isspace(c); }).base();
  if (first >= last)
    return "";
  return std::string(first, last);
}

std::string stripComment(const std::string &line) {
  const auto pos = line.find('#');
  return trim(pos == std::string::npos ? line : line.substr(0, pos));
}

std::string csvEscape(const std::string &text) {
  if (text.find_first_of(",\"\n") == std::string::npos)
    return text;
  std::string out = "\"";
  for (char c : text) {
    if (c == '"')
      out += "\"\"";
    else
      out += c;
  }
  out += "\"";
  return out;
}

std::string xmlEscape(const std::string &text) {
  std::string out;
  for (char c : text) {
    switch (c) {
      case '&': out += "&amp;"; break;
      case '"': out += "&quot;"; break;
      case '<': out += "&lt;"; break;
      case '>': out += "&gt;"; break;
      default: out += c; break;
    }
  }
  return out;
}

std::string shellQuote(const std::string &text) {
  std::string out = "'";
  for (char c : text) {
    if (c == '\'')
      out += "'\\''";
    else
      out += c;
  }
  out += "'";
  return out;
}

std::string joinModules(const std::vector<std::string> &modules) {
  std::ostringstream out;
  for (size_t i = 0; i < modules.size(); i++) {
    if (i)
      out << ' ';
    out << modules[i];
  }
  return out.str();
}

std::string zeroList(int count) {
  std::ostringstream out;
  for (int i = 0; i < count; i++) {
    if (i)
      out << ", ";
    out << '0';
  }
  return out.str();
}

std::string sanitizeName(std::string name) {
  for (char &c : name) {
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-')
      c = '_';
  }
  return name;
}

bool hasPathPart(const fs::path &path, const std::string &part) {
  for (const auto &pathPart : path) {
    if (pathPart == part)
      return true;
  }
  return false;
}

bool isUnsupportedStandaloneUrdf(const fs::path &urdfPath) {
  const std::string filename = urdfPath.filename().string();
  if (filename == "trackedTemplate.urdf")
    return true;
  if (hasPathPart(urdfPath, "raw"))
    return true;
  if (filename == "raipal_base.urdf")
    return true;
  if (filename == "raipal_upper-only_L.urdf" || filename == "raipal_upper-only_R.urdf")
    return true;
  return false;
}

bool isUnsupportedModuleCase(const ModuleCase &moduleCase) {
  return moduleCase.target == "raipal_base" ||
         moduleCase.target == "raipal_upper-only_L" ||
         moduleCase.target == "raipal_upper-only_R";
}

template <typename Derived>
std::string eigenValues(const Eigen::MatrixBase<Derived> &values, int maxItems = 96) {
  std::ostringstream out;
  out << std::setprecision(12);
  const Eigen::Index n = values.size();
  const auto eval = values.derived().eval();
  for (Eigen::Index i = 0; i < n && i < maxItems; i++) {
    if (i)
      out << ' ';
    out << eval(i);
  }
  if (n > maxItems)
    out << " ...";
  return out.str();
}

std::string stringValues(const std::vector<std::string> &values) {
  std::ostringstream out;
  for (size_t i = 0; i < values.size(); i++) {
    if (i)
      out << " | ";
    out << i << ':' << values[i];
  }
  return out.str();
}

class Reporter {
 public:
  Reporter(const fs::path &csvPath, const fs::path &txtPath, bool writeHeader)
      : csv_(csvPath), txt_(txtPath) {
    if (!csv_)
      throw std::runtime_error("failed to open CSV output: " + csvPath.string());
    if (!txt_)
      throw std::runtime_error("failed to open TXT output: " + txtPath.string());
    if (writeHeader)
      csv_ << "scenario,source,status,metric,size,expected,delta,values,notes\n";
  }

  template <typename Derived>
  void row(const std::string &scenario,
           const std::string &source,
           const std::string &metric,
           Eigen::Index size,
           Eigen::Index expected,
           const Eigen::MatrixBase<Derived> &values,
           const std::string &notes = "") {
    rowText(scenario, source, "ok", metric, size, expected, eigenValues(values), notes);
  }

  void rowSize(const std::string &scenario,
               const std::string &source,
               const std::string &metric,
               long long size,
               long long expected,
               const std::string &values = "",
               const std::string &notes = "") {
    rowText(scenario, source, "ok", metric, size, expected, values, notes);
  }

  void error(const std::string &scenario, const std::string &source, const std::string &message) {
    rowText(scenario, source, "error", "load", -1, -1, "", message);
  }

  std::ostream &txt() { return txt_; }

 private:
  void rowText(const std::string &scenario,
               const std::string &source,
               const std::string &status,
               const std::string &metric,
               long long size,
               long long expected,
               const std::string &values,
               const std::string &notes) {
    const long long delta = (size >= 0 && expected >= 0) ? size - expected : 0;
    csv_ << csvEscape(scenario) << ','
         << csvEscape(source) << ','
         << csvEscape(status) << ','
         << csvEscape(metric) << ','
         << size << ','
         << expected << ','
         << delta << ','
         << csvEscape(values) << ','
         << csvEscape(notes) << '\n';

    txt_ << "[" << status << "] " << scenario << "\n"
         << "  source: " << source << "\n"
         << "  " << metric << ": size=" << size << ", expected=" << expected
         << ", delta=" << delta << "\n";
    if (!values.empty())
      txt_ << "  values: " << values << "\n";
    if (!notes.empty())
      txt_ << "  notes: " << notes << "\n";
    txt_ << "\n";
  }

  std::ofstream csv_;
  std::ofstream txt_;
};

void reportSystem(Reporter &reporter,
                  const std::string &scenario,
                  const std::string &source,
                  raisim::ArticulatedSystem *robot) {
  if (!robot) {
    reporter.error(scenario, source, "robot pointer was null");
    return;
  }

  const auto dof = static_cast<long long>(robot->getDOF());
  const auto gcDim = static_cast<long long>(robot->getGeneralizedCoordinateDim());
  const auto gvDim = static_cast<long long>(robot->getGeneralizedVelocityDim());

  reporter.txt() << "========== " << scenario << " ==========\n"
                 << "source: " << source << "\n"
                 << "DOF=" << dof << ", gcDim=" << gcDim << ", gvDim=" << gvDim << "\n\n";

  reporter.rowSize(scenario, source, "getDOF", dof, dof);
  reporter.rowSize(scenario, source, "getGeneralizedCoordinateDim", gcDim, gcDim);
  reporter.rowSize(scenario, source, "getGeneralizedVelocityDim", gvDim, gvDim);

  reporter.row(scenario, source, "getGeneralizedCoordinate",
               robot->getGeneralizedCoordinate().e().size(), gcDim, robot->getGeneralizedCoordinate().e());
  reporter.row(scenario, source, "getGeneralizedVelocity",
               robot->getGeneralizedVelocity().e().size(), gvDim, robot->getGeneralizedVelocity().e());
  reporter.row(scenario, source, "getGeneralizedAcceleration",
               robot->getGeneralizedAcceleration().e().size(), gvDim, robot->getGeneralizedAcceleration().e());
  reporter.row(scenario, source, "getGeneralizedForce",
               robot->getGeneralizedForce().e().size(), gvDim, robot->getGeneralizedForce().e());
  reporter.row(scenario, source, "getFeedForwardGeneralizedForce",
               robot->getFeedForwardGeneralizedForce().e().size(), gvDim, robot->getFeedForwardGeneralizedForce().e());

  Eigen::VectorXd gc;
  Eigen::VectorXd gv;
  robot->getState(gc, gv);
  reporter.row(scenario, source, "getState.genco", gc.size(), gcDim, gc);
  reporter.row(scenario, source, "getState.genvel", gv.size(), gvDim, gv);

  Eigen::VectorXd pTarget = Eigen::VectorXd::Zero(gcDim);
  Eigen::VectorXd dTarget = Eigen::VectorXd::Zero(gvDim);
  robot->getPdTarget(pTarget, dTarget);
  reporter.row(scenario, source, "getPdTarget.position", pTarget.size(), gcDim, pTarget);
  reporter.row(scenario, source, "getPdTarget.velocity", dTarget.size(), gvDim, dTarget);

  Eigen::VectorXd pGain;
  Eigen::VectorXd dGain;
  robot->getPdGains(pGain, dGain);
  reporter.row(scenario, source, "getPdGains.p", pGain.size(), gvDim, pGain);
  reporter.row(scenario, source, "getPdGains.d", dGain.size(), gvDim, dGain);

  reporter.row(scenario, source, "getActuationUpperLimits",
               robot->getActuationUpperLimits().e().size(), gvDim, robot->getActuationUpperLimits().e());
  reporter.row(scenario, source, "getActuationLowerLimits",
               robot->getActuationLowerLimits().e().size(), gvDim, robot->getActuationLowerLimits().e());
  reporter.row(scenario, source, "getRotorInertia",
               robot->getRotorInertia().e().size(), gvDim, robot->getRotorInertia().e());
  reporter.row(scenario, source, "getJointVelocityLimits",
               robot->getJointVelocityLimits().e().size(), gvDim, robot->getJointVelocityLimits().e(),
               "expected column is gvDim; nonzero delta means this getter is not gv-indexed");

  const auto &jointLimits = robot->getJointLimits();
  reporter.rowSize(scenario, source, "getJointLimits",
                   static_cast<long long>(jointLimits.size()), -1, "", "std::vector<Vec<2>>; joint-list indexed");
  reporter.rowSize(scenario, source, "getMovableJointNames",
                   static_cast<long long>(robot->getMovableJointNames().size()),
                   static_cast<long long>(robot->getMovableJointNames().size()),
                   stringValues(robot->getMovableJointNames()));
  reporter.rowSize(scenario, source, "getFrames", static_cast<long long>(robot->getFrames().size()), -1);
  reporter.rowSize(scenario, source, "getMass", static_cast<long long>(robot->getMass().size()), -1);
  reporter.rowSize(scenario, source, "getInertia", static_cast<long long>(robot->getInertia().size()), -1);
  reporter.rowSize(scenario, source, "getBodyCOM_B", static_cast<long long>(robot->getBodyCOM_B().size()), -1);
  reporter.rowSize(scenario, source, "getBodyCOM_W", static_cast<long long>(robot->getBodyCOM_W().size()), -1);
  reporter.rowSize(scenario, source, "getJointPos_P", static_cast<long long>(robot->getJointPos_P().size()), -1);
  reporter.rowSize(scenario, source, "getJointOrientation_P", static_cast<long long>(robot->getJointOrientation_P().size()), -1);
  reporter.rowSize(scenario, source, "getJointAxis_P", static_cast<long long>(robot->getJointAxis_P().size()), -1);

  std::ostringstream map;
  const auto &names = robot->getMovableJointNames();
  for (size_t i = 0; i < names.size(); i++) {
    if (i)
      map << " | ";
    map << i << ':' << names[i] << "->gv" << robot->getGeneralizedVelocityIndex(names[i]);
  }
  reporter.rowSize(scenario, source, "movableJointNameToGvIndex",
                   static_cast<long long>(names.size()), static_cast<long long>(names.size()), map.str());
}

fs::path resourceRootFromBuild() {
#ifdef RESOURCE_DIR
  fs::path root(STRINGIFY(RESOURCE_DIR));
  if (fs::exists(root))
    return fs::canonical(root);
#endif
  fs::path cwdRoot = fs::current_path() / "../resource";
  if (fs::exists(cwdRoot))
    return fs::canonical(cwdRoot);
  return fs::canonical(fs::current_path() / "resource");
}

std::vector<fs::path> findUrdfs(const fs::path &resourceRoot) {
  std::vector<fs::path> urdfs;
  for (const auto &entry : fs::recursive_directory_iterator(resourceRoot)) {
    if (!entry.is_regular_file())
      continue;
    bool generatedByThisChecker = false;
    for (const auto &part : entry.path()) {
      if (part == "parameter_length_check" || part == "parameter_length_check_generated_bases") {
        generatedByThisChecker = true;
        break;
      }
    }
    if (generatedByThisChecker)
      continue;
    if (entry.path().extension() != ".urdf")
      continue;
    if (entry.path().filename().string().rfind("generated_", 0) == 0)
      continue;
    if (isUnsupportedStandaloneUrdf(entry.path()))
      continue;
    urdfs.push_back(fs::canonical(entry.path()));
  }
  std::sort(urdfs.begin(), urdfs.end());
  return urdfs;
}

std::vector<fs::path> findSeedFiles(const fs::path &resourceRoot) {
  std::vector<fs::path> seeds;
  for (const auto &entry : fs::recursive_directory_iterator(resourceRoot)) {
    if (!entry.is_regular_file())
      continue;
    bool generatedByThisChecker = false;
    for (const auto &part : entry.path()) {
      if (part == "parameter_length_check" || part == "parameter_length_check_generated_bases") {
        generatedByThisChecker = true;
        break;
      }
    }
    if (generatedByThisChecker)
      continue;
    const auto filename = entry.path().filename().string();
    const auto extension = entry.path().extension().string();
    if (filename.find("seed") != std::string::npos && (extension == ".yaml" || extension == ".yml" || extension == ".urdf"))
      seeds.push_back(fs::canonical(entry.path()));
  }
  std::sort(seeds.begin(), seeds.end());
  return seeds;
}

std::vector<ModuleCase> parseSeedFile(const fs::path &seedPath) {
  std::ifstream input(seedPath);
  if (!input)
    throw std::runtime_error("failed to open seed file: " + seedPath.string());

  std::vector<ModuleCase> cases;
  ModuleCase current;
  current.seedPath = seedPath;
  current.packageRoot = seedPath.parent_path().parent_path();

  std::string line;
  while (std::getline(input, line)) {
    const std::string stripped = stripComment(line);
    if (stripped.empty())
      continue;

    if (stripped.back() == ':' && stripped.front() != '-') {
      if (!current.target.empty() && !current.modules.empty())
        cases.push_back(current);
      current = ModuleCase{};
      current.seedPath = seedPath;
      current.packageRoot = seedPath.parent_path().parent_path();
      current.target = trim(stripped.substr(0, stripped.size() - 1));
      continue;
    }

    if (!stripped.empty() && stripped.front() == '-') {
      std::string module = trim(stripped.substr(1));
      if (!module.empty())
        current.modules.push_back(module);
    }
  }

  if (!current.target.empty() && !current.modules.empty())
    cases.push_back(current);
  return cases;
}

std::vector<ModuleCase> findModuleCases(const fs::path &resourceRoot) {
  std::vector<ModuleCase> cases;
  std::set<std::string> seen;
  for (const auto &seed : findSeedFiles(resourceRoot)) {
    for (const auto &moduleCase : parseSeedFile(seed)) {
      if (isUnsupportedModuleCase(moduleCase))
        continue;
      std::ostringstream key;
      key << moduleCase.packageRoot << '|' << moduleCase.target << '|' << joinModules(moduleCase.modules);
      if (seen.insert(key.str()).second)
        cases.push_back(moduleCase);
    }
  }
  std::sort(cases.begin(), cases.end(), [](const ModuleCase &a, const ModuleCase &b) {
    if (a.packageRoot != b.packageRoot)
      return a.packageRoot < b.packageRoot;
    return a.target < b.target;
  });
  return cases;
}

fs::path findOrCreateModuleBaseUrdf(const ModuleCase &moduleCase,
                                    const fs::path &generatedBaseRoot,
                                    fs::path *urdfPathForXml) {
  const std::string package = moduleCase.packageRoot.filename().string();
  std::vector<fs::path> candidates = {
      moduleCase.packageRoot / "urdf" / (package + "_base.urdf"),
      moduleCase.packageRoot / "urdf" / "raipal_base.urdf",
      moduleCase.packageRoot / "urdf" / "head_base.urdf",
      moduleCase.packageRoot / "urdf" / (moduleCase.target + "_base.urdf"),
      moduleCase.packageRoot / "urdf" / "base.urdf",
  };
  const fs::path urdfDir = moduleCase.packageRoot / "urdf";
  if (fs::exists(urdfDir)) {
    for (const auto &entry : fs::directory_iterator(urdfDir)) {
      if (entry.is_regular_file() && entry.path().filename().string().find("_base.urdf") != std::string::npos)
        candidates.push_back(entry.path());
    }
  }

  for (const auto &candidate : candidates) {
    if (fs::exists(candidate)) {
      *urdfPathForXml = fs::canonical(candidate);
      return candidate;
    }
  }

  if (moduleCase.modules.empty())
    throw std::runtime_error("module case has no modules: " + moduleCase.target);

  const fs::path modulePath = moduleCase.packageRoot / "modules" / moduleCase.modules.front();
  if (!fs::exists(modulePath))
    throw std::runtime_error("first module file is missing: " + modulePath.string());

  const fs::path generatedDir = generatedBaseRoot / sanitizeName(package);
  fs::create_directories(generatedDir);
  const fs::path generatedBase = generatedDir / (sanitizeName(moduleCase.target + "_" + moduleCase.modules.front()) + ".urdf");

  std::ifstream module(modulePath);
  std::ofstream urdf(generatedBase);
  if (!module)
    throw std::runtime_error("failed to open module file: " + modulePath.string());
  if (!urdf)
    throw std::runtime_error("failed to write generated base URDF: " + generatedBase.string());

  urdf << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
       << "<robot name=\"" << xmlEscape(moduleCase.target + "_base") << "\">\n"
       << module.rdbuf()
       << "\n</robot>\n";

  *urdfPathForXml = fs::canonical(generatedBase);
  return generatedBase;
}

GeneratedModuleConfig writeModuleConfig(const ModuleCase &moduleCase,
                                        const fs::path &configDir,
                                        const fs::path &generatedBaseRoot) {
  fs::create_directories(configDir);
  const std::string package = moduleCase.packageRoot.filename().string();
  const fs::path configPath = configDir / (sanitizeName(package + "_" + moduleCase.target) + ".xml");
  fs::path urdfPathForXml;
  findOrCreateModuleBaseUrdf(moduleCase, generatedBaseRoot, &urdfPathForXml);

  std::vector<std::string> remainingModules = moduleCase.modules;
  if (!remainingModules.empty() && !remainingModules.front().empty() && remainingModules.front()[0] == '_')
    remainingModules.erase(remainingModules.begin());

  std::ofstream xml(configPath);
  if (!xml)
    throw std::runtime_error("failed to write generated config: " + configPath.string());

  xml << "<?xml version=\"1.0\" ?>\n"
      << "<raisim version=\"1.0\">\n"
      << "  <gravity value=\"0.000000, 0.000000, -9.810000\" />\n"
      << "  <timeStep value=\"0.001\" />\n"
      << "  <objects>\n"
      << "    <articulatedSystem name=\"robot\"\n"
      << "      resDir=\"" << xmlEscape(moduleCase.packageRoot.string()) << "\"\n"
      << "      urdfPath=\"" << xmlEscape(urdfPathForXml.string()) << "\"\n"
      << "      collisionGroup=\"1\" collisionMask=\"-1\"\n"
      << "      modules=\"@@MODULES\">\n"
      << "      <state qpos=\"" << zeroList(64) << "\" qvel=\"" << zeroList(64) << "\" />\n"
      << "    </articulatedSystem>\n"
      << "  </objects>\n"
      << "</raisim>\n";
  return {configPath, joinModules(remainingModules)};
}

void initializeOutputs(const fs::path &csvPath, const fs::path &txtPath) {
  std::ofstream csv(csvPath);
  std::ofstream txt(txtPath);
  if (!csv)
    throw std::runtime_error("failed to open CSV output: " + csvPath.string());
  if (!txt)
    throw std::runtime_error("failed to open TXT output: " + txtPath.string());
  csv << "scenario,source,status,metric,size,expected,delta,values,notes\n";
}

void appendChildOutput(const fs::path &finalCsv,
                       const fs::path &finalTxt,
                       const fs::path &childCsv,
                       const fs::path &childTxt,
                       const fs::path &childLog) {
  {
    std::ifstream in(childCsv);
    std::ofstream out(finalCsv, std::ios::app);
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
      if (first) {
        first = false;
        continue;
      }
      out << line << '\n';
    }
  }
  {
    std::ifstream in(childTxt);
    std::ofstream out(finalTxt, std::ios::app);
    out << in.rdbuf();
    std::ifstream log(childLog);
    if (log.peek() != std::ifstream::traits_type::eof()) {
      out << "----- child stderr/stdout: " << childLog << " -----\n";
      out << log.rdbuf() << "\n";
    }
  }
}

void appendErrorRow(const fs::path &finalCsv,
                    const fs::path &finalTxt,
                    const std::string &scenario,
                    const std::string &source,
                    const std::string &message) {
  std::ofstream csv(finalCsv, std::ios::app);
  std::ofstream txt(finalTxt, std::ios::app);
  csv << csvEscape(scenario) << ','
      << csvEscape(source) << ",error,load,-1,-1,0,,"
      << csvEscape(message) << '\n';
  txt << "[error] " << scenario << "\n"
      << "  source: " << source << "\n"
      << "  notes: " << message << "\n\n";
}

int runChild(const fs::path &exe,
             const std::string &mode,
             const std::string &scenario,
             const fs::path &source,
             const fs::path &childPrefix) {
  const fs::path logPath = childPrefix.string() + ".log";
  const std::string command =
      shellQuote(exe.string()) + " " +
      shellQuote(mode) + " " +
      shellQuote(scenario) + " " +
      shellQuote(source.string()) + " " +
      shellQuote(childPrefix.string()) + " > " +
      shellQuote(logPath.string()) + " 2>&1";
  return std::system(command.c_str());
}

int runModuleChild(const fs::path &exe,
                   const std::string &scenario,
                   const fs::path &source,
                   const std::string &modules,
                   const fs::path &childPrefix) {
  const fs::path logPath = childPrefix.string() + ".log";
  const std::string command =
      shellQuote(exe.string()) + " " +
      shellQuote("--single-module-config") + " " +
      shellQuote(scenario) + " " +
      shellQuote(source.string()) + " " +
      shellQuote(modules) + " " +
      shellQuote(childPrefix.string()) + " > " +
      shellQuote(logPath.string()) + " 2>&1";
  return std::system(command.c_str());
}

void runSingleDirect(const std::string &scenario, const fs::path &urdfPath, const fs::path &outputPrefix) {
  Reporter reporter(outputPrefix.string() + ".csv", outputPrefix.string() + ".txt", true);
  try {
    raisim::World world;
    auto *robot = world.addArticulatedSystem(urdfPath.string());
    reportSystem(reporter, scenario, urdfPath.string(), robot);
  } catch (const std::exception &e) {
    reporter.error(scenario, urdfPath.string(), e.what());
    throw;
  }
}

void runSingleConfig(const std::string &scenario, const fs::path &configPath, const fs::path &outputPrefix) {
  Reporter reporter(outputPrefix.string() + ".csv", outputPrefix.string() + ".txt", true);
  try {
    raisim::World world(configPath.string());
    auto *robot = reinterpret_cast<raisim::ArticulatedSystem*>(world.getObject("robot"));
    reportSystem(reporter, scenario, configPath.string(), robot);
  } catch (const std::exception &e) {
    reporter.error(scenario, configPath.string(), e.what());
    throw;
  }
}

void runSingleModuleConfig(const std::string &scenario,
                           const fs::path &configPath,
                           const std::string &modules,
                           const fs::path &outputPrefix) {
  Reporter reporter(outputPrefix.string() + ".csv", outputPrefix.string() + ".txt", true);
  try {
    raisim::World world(configPath.string(), {{"MODULES", modules}});
    auto *robot = reinterpret_cast<raisim::ArticulatedSystem*>(world.getObject("robot"));
    reportSystem(reporter, scenario, configPath.string() + " MODULES=" + modules, robot);
  } catch (const std::exception &e) {
    reporter.error(scenario, configPath.string() + " MODULES=" + modules, e.what());
    throw;
  }
}

int runAll(const fs::path &exe, const fs::path &outputRoot) {
  const fs::path resourceRoot = resourceRootFromBuild();
  const fs::path sandboxRoot = resourceRoot.parent_path();
  const fs::path finalCsv = outputRoot / "parameter_length_check.csv";
  const fs::path finalTxt = outputRoot / "parameter_length_check.txt";
  const fs::path childDir = outputRoot / "children-logs";
  const fs::path configDir = outputRoot / "generated-configs";
  const fs::path generatedBaseDir = outputRoot / "generated-bases";

  fs::create_directories(outputRoot);
  if (fs::exists(childDir))
    fs::remove_all(childDir);
  if (fs::exists(configDir))
    fs::remove_all(configDir);
  if (fs::exists(generatedBaseDir))
    fs::remove_all(generatedBaseDir);
  fs::create_directories(childDir);
  fs::create_directories(configDir);
  fs::create_directories(generatedBaseDir);
  initializeOutputs(finalCsv, finalTxt);

  {
    std::ofstream txt(finalTxt, std::ios::app);
    txt << "resourceRoot: " << resourceRoot << "\n"
        << "sandboxRoot: " << sandboxRoot << "\n"
        << "outputRoot: " << outputRoot << "\n"
        << "childrenLogs: " << childDir << "\n"
        << "generatedConfigs: " << configDir << "\n"
        << "generatedBases: " << generatedBaseDir << "\n\n";
  }

  const auto urdfs = findUrdfs(resourceRoot);
  const auto moduleCases = findModuleCases(resourceRoot);
  {
    std::ofstream txt(finalTxt, std::ios::app);
    txt << "discovered direct URDFs: " << urdfs.size() << "\n"
        << "discovered module seed cases: " << moduleCases.size() << "\n\n";
  }

  size_t childIndex = 0;
  int failures = 0;

  for (const auto &urdf : urdfs) {
    const auto rel = fs::relative(urdf, resourceRoot).string();
    const std::string scenario = "direct urdf: " + rel;
    const fs::path childPrefix = childDir / ("direct_" + std::to_string(childIndex++));
    const int status = runChild(exe, "--single-direct", scenario, urdf, childPrefix);
    appendChildOutput(finalCsv, finalTxt, childPrefix.string() + ".csv", childPrefix.string() + ".txt", childPrefix.string() + ".log");
    if (status != 0) {
      failures++;
      appendErrorRow(finalCsv, finalTxt, scenario, urdf.string(), "child process returned status " + std::to_string(status));
    }
  }

  for (const auto &moduleCase : moduleCases) {
    const auto relRoot = fs::relative(moduleCase.packageRoot, resourceRoot).string();
    const std::string scenario = "module config: " + relRoot + "/" + moduleCase.target;
    const auto generated = writeModuleConfig(moduleCase, configDir, generatedBaseDir);
    {
      std::ofstream txt(finalTxt, std::ios::app);
      txt << "generated config: " << generated.configPath << "\n"
          << "  seed: " << moduleCase.seedPath << "\n"
          << "  modules from seed: " << joinModules(moduleCase.modules) << "\n"
          << "  modules passed to World: " << generated.modulesForWorld << "\n\n";
    }
    const fs::path childPrefix = childDir / ("module_" + std::to_string(childIndex++));
    const int status = runModuleChild(exe, scenario, generated.configPath, generated.modulesForWorld, childPrefix);
    appendChildOutput(finalCsv, finalTxt, childPrefix.string() + ".csv", childPrefix.string() + ".txt", childPrefix.string() + ".log");
    if (status != 0) {
      failures++;
      appendErrorRow(finalCsv, finalTxt, scenario, generated.configPath.string(), "child process returned status " + std::to_string(status));
    }
  }

  std::cout << "wrote " << finalCsv << "\n"
            << "wrote " << finalTxt << "\n"
            << "wrote child logs under " << childDir << "\n"
            << "wrote generated configs under " << configDir << "\n"
            << "wrote generated bases under " << generatedBaseDir << "\n"
            << "direct URDF cases: " << urdfs.size() << "\n"
            << "module config cases: " << moduleCases.size() << "\n"
            << "child failures: " << failures << "\n";

  return failures == 0 ? 0 : 1;
}

}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc == 5 && std::string(argv[1]) == "--single-direct") {
      runSingleDirect(argv[2], argv[3], argv[4]);
      return 0;
    }
    if (argc == 5 && std::string(argv[1]) == "--single-config") {
      runSingleConfig(argv[2], argv[3], argv[4]);
      return 0;
    }
    if (argc == 6 && std::string(argv[1]) == "--single-module-config") {
      runSingleModuleConfig(argv[2], argv[3], argv[4], argv[5]);
      return 0;
    }

    const fs::path resourceRoot = resourceRootFromBuild();
    const fs::path defaultOutputRoot = resourceRoot / "parameter_length_check";
    const fs::path outputRoot = argc > 1 ? fs::path(argv[1]) : defaultOutputRoot;
    return runAll(fs::canonical(argv[0]), outputRoot);
  } catch (const std::exception &e) {
    std::cerr << "parameter_length_check failed: " << e.what() << "\n";
    return 2;
  }
}
